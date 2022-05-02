void compositeEvent() 
{
  int station = 21; 

  int deep_run = 345; 
  int deep_evt = 0; 
  mattak::Dataset ds_deep(station,deep_run,0); 
  ds_deep.setEntry(deep_evt); 
  mattak::WaveformPlotOptions deep_opt; 
  deep_opt.mask = 0xe00f1f; 
  deep_opt.autoscale = false;
  deep_opt.min=-900;
  deep_opt.max=900;

  int shallow_run = 357; 
  int shallow_evt = 0; 
  mattak::Dataset ds_shallow(station,shallow_run,0); 
  ds_shallow.setEntry(shallow_evt); 
  mattak::WaveformPlotOptions shallow_opt; 
  shallow_opt.mask = 0x1ff000; 
  shallow_opt.autoscale = false;
  shallow_opt.min=-900;
  shallow_opt.max=900;


  int intermed_run = 351; 
  int intermed_evt = 0;
  mattak::Dataset ds_intermed(station,intermed_run,0); 
  ds_intermed.setEntry(intermed_evt); 
  mattak::WaveformPlotOptions intermed_opt; 
  intermed_opt.mask = 0xe0; 
  intermed_opt.autoscale = false;
  intermed_opt.min=-900;
  intermed_opt.max=900;



  TCanvas * call = new TCanvas("call","call",2400,1300); 
  call->cd() ;
  TPad * pdeep = new TPad("pdeep","deep", 0,0.0,1,0.5,kWhite);
  pdeep->Draw();
  ds_deep.raw()->drawWaveforms(deep_opt, pdeep);
  call->cd() ;
  TPad * pshallow = new TPad("pshallow","shallow", 0,0.625,1,1);
  pshallow->Draw();
  ds_shallow.raw()->drawWaveforms(shallow_opt, pshallow);
  call->cd() ;
  TPad * pintermed = new TPad("pintermed","intermed", 0.,0.5,1,0.625,kWhite); 
  pintermed->Draw(); 
  ds_intermed.raw()->drawWaveforms(intermed_opt, pintermed);



}
