
#include "mattak/Detector.h"
#include <iostream>
#include <map>

#ifdef MONGO_SUPPORT

#include <mongoc/mongoc.h>
#include <bson/bson.h>
#include "TMutex.h" 

static bool mongo_init;
static TMutex mongo_lock;

//RAII wrapper
struct MongoClient
{
  std::string connection_string;
  mongoc_client_t * client  = nullptr;
  mongoc_database_t * db  = nullptr;
  bool ok = false;

  MongoClient (const char * conn, const char * dbname = "RNOG_Live")
    : connection_string(conn)
  {

    if (!mongo_init)
    {
      TLockGuard l(&mongo_lock);
      if (!mongo_init)
      {
        mongoc_init();
        mongo_init = true;
      }
    }

    client = mongoc_client_new(conn);
    if (!client)
    {
      std::cerr<< "Could not open " << conn << std::endl;
      return;
    }

    db = mongoc_client_get_database(client, dbname);

     if (!db)
    {
      std::cerr<< "Could not find " << dbname << std::endl;
      return;
    }
     ok = true;
  }


  ~MongoClient()
  {
    if (db) mongoc_database_destroy(db);
    if (client) mongoc_client_destroy(client);
  }

};

std::map<const char *,MongoClient> clients;



#define BSON_IT_GET(X,type,it)  if (!strcmp(key,#X)) { X = bson_iter_##type(it); }
#define BSON_IT_GET_STR(X,it)  if (!strcmp(key,#X)) { X = bson_iter_utf8(it,0); }
#define BSON_IT_GET_TIME(X,it)  if (!strcmp(key,#X)) { int64_t t  =bson_iter_date_time(it); X = TTimeStamp(t/1000, 1e6 * (t % 1000)); }

bool setTVector(bson_iter_t *it, TVector3 * v)
{
  bson_iter_t a;
  bson_iter_recurse(it,&a);
  int ipos = 0;

  while (bson_iter_next(&a))
  {
    if (ipos > 2)
    {
      std::cerr << "Expected only 3 elements " << std::endl;
      return false;
    }

    v->operator[](ipos++) = bson_iter_as_double(&a);
  }

  if (ipos < 3)
  {
    std::cerr << "Expected 3 elements" << std::endl;
    return false;
  }

  return true;
}

struct DBChannelInfo
{
  DBChannelInfo(bson_iter_t *it)
  {
    while (bson_iter_next(it))
    {
      const char * key = bson_iter_key(it);
      std::string ant_type;
      BSON_IT_GET(id,int32, it)
      BSON_IT_GET_TIME(commission_time, it)
      BSON_IT_GET_TIME(decommission_time, it)
      BSON_IT_GET_STR(id_position, it)
      BSON_IT_GET_STR(id_signal, it)
      BSON_IT_GET_STR(ant_type, it);
      type = toupper(ant_type[0]);
      if (!strcmp(key,"installed_components"))
      {
        bson_iter_t ic_it;
        bson_iter_recurse(it,&ic_it);
        while (bson_iter_next(&ic_it)) installed_components[bson_iter_key(&ic_it)]=bson_iter_utf8(&ic_it,0);
      }
    }
  }

  int id = -1;
  TTimeStamp commission_time;
  TTimeStamp decommission_time;
  std::string id_position;
  std::string id_signal;
  char type='\0'; //H,V,L
  std::map<std::string,std::string> installed_components;
};


const mattak::Detector* mattak::Detector::fromDB(int station, const TTimeStamp & t, const char * connection_string)
{
  if (!clients.count(connection_string))
  {
    TLockGuard lck(&mongo_lock);

    if (!clients.count(connection_string))
    {
      clients.emplace(connection_string, connection_string);
    }
  }


  MongoClient & client = (*clients.find(connection_string)).second;
  if (!client.ok)
  {
    std::cerr << "MongoDB client not ok for " << connection_string << std::endl;
    return nullptr;
  }

  // get general station information

  mongoc_collection_t * station_rnog = mongoc_database_get_collection(client.db,"station_rnog");

  if (!station_rnog)
  {
    std::cerr << "Cannot find station_rnog in database" << std::endl;
    return nullptr;
  }

  double now = t.AsDouble() *1000;
  // just cargo-culted from NMC
  bson_t * time_filter = BCON_NEW("pipeline","["
                                    "$match",
                                    "{",
                                       "commission_time", "{", "$lte", BCON_DATE_TIME(now), "}",
                                       "commission_time", "{", "$gte", BCON_DATE_TIME(now), "}",
                                       "id", BCON_INT32(station),
                                     "}",
                                    "]");


  mongoc_cursor_t *c =  mongoc_collection_aggregate(station_rnog, MONGOC_QUERY_NONE, time_filter,  nullptr, nullptr);

  bson_destroy(time_filter);
  const bson_t * result;
  int nresults = 0;


  //station position id
  std::string id_position = "";

  std::vector<DBChannelInfo> channel_infos;;
  while (mongoc_cursor_next(c, &result))
  {

    if (nresults++)
    {
      std::cerr << "WARNING: MULTIPLE RESULTS WERE FOUND FOR A GIVEN TIME. I AM JUST USING THE FIRST ONE" << std::endl;
      break;
    }

    bson_iter_t it;
    bson_iter_init(&it, result);
    while (bson_iter_next(&it))
    {
      const char * key = bson_iter_key(&it);
      BSON_IT_GET_STR(id_position,&it)

      if (!strcmp(key,"channels"))
      {
        bson_iter_t ch_it;
        bson_iter_recurse(&it,&ch_it);

        channel_infos.emplace_back(&ch_it);

        //reject if out of bounds
        if (t.AsDouble() > channel_infos.back().decommission_time.AsDouble() ||  t.AsDouble() < channel_infos.back().commission_time.AsDouble())
        {
          channel_infos.pop_back();
        }
      }
    }
  }

  //cleanup station_rnog
  mongoc_cursor_destroy(c);
  mongoc_collection_destroy(station_rnog);

  if (!nresults)
  {
    std::cerr << "No detector found for station/time" << std::endl;
    return nullptr;
  }


  Detector * d = new Detector();
  d->det_time = t;

  //get station_position
  mongoc_collection_t * station_position = mongoc_database_get_collection(client.db,"station_position");
  if (!station_position)
  {
    std::cerr << "Can't find station position" << std::endl;
    delete d;
    return nullptr;
  }
  bson_t* station_position_query = BCON_NEW("$query","{","id",BCON_UTF8(id_position.c_str()),"}");

  c = mongoc_collection_find_with_opts(station_position, station_position_query,NULL,NULL);
  bson_destroy(station_position_query);

  nresults = 0;
  while (mongoc_cursor_next(c, &result))
  {
    if (nresults++)
    {
      std::cerr << "WARNING: MULTIPLE RESULTS FOUND FOR " << id_position << std::endl;
      break;
    }
    bson_iter_t it;
    bson_iter_init(&it,result);
    bson_iter_find(&it,"position");
    setTVector(&it, &d->station_position);
  };


  mongoc_cursor_destroy(c);
  mongoc_collection_destroy(station_position);



  //query channel positions
  mongoc_collection_t * channel_position = mongoc_database_get_collection(client.db, "channel_position");
  bson_t * channel_query = bson_new();
  bson_t * channel_info_ids;
  bson_append_document_begin(channel_query,"id",2,channel_info_ids);

  bson_array_builder_t *ab;
  bson_append_array_builder_begin(channel_info_ids,"$in",3,&ab);
  //build the query
  for (auto c : channel_infos)
  {
    bson_array_builder_append_utf8(ab, c.id_position.c_str(), c.id_position.size());
  }
  bson_append_array_builder_end(channel_info_ids,ab);
  bson_append_document_end(channel_query,channel_info_ids);


  c = mongoc_collection_find_with_opts(channel_position,channel_query,NULL,NULL);
  if (!c) 
  {
    std::cerr << "Could not find channel positions" << std::endl;
    delete d;
    return nullptr;

  }
  bson_destroy(channel_query);

  nresults = 0;
  while(mongoc_cursor_next(c,&result))
  {
    bson_iter_t it;
    bson_iter_init(&it,result);

    int channel_id;
    TVector3 v;
    while(bson_iter_next(&it))
    {
      const char * key = bson_iter_key(&it);

      BSON_IT_GET(channel_id,int32, &it);
      if (!strcmp(key,"position")) setTVector(&it,&v);
    }
    if (channel_id > mattak::k::num_radiant_channels)
    {
      std::cerr << "Got channel id bigger than 23..." << std::endl;
      continue;
    }

    d->antenna_positions[channel_id] = v;
    d->global_antenna_positions[channel_id] = v + d->station_position;
  }

  mongoc_cursor_destroy(c);
  mongoc_collection_destroy(channel_position);

  return d;
}

#else

const mattak::Detector* mattak::Detector::fromDB(int station, const TTimeStamp & t, const char * connection_string)
{
  (void) station;
  (void) t;
  (void) connection_string;
  std::cerr << "mattak compiled without mongodb support " << std::endl;
  return nullptr;
}


#endif

