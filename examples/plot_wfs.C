
void plot_wfs(const char * infile, int i = 0, const char * tree_name = "waveforms", const char * branch_name = "waveforms")
{
  TFile f(infile); 
  TTree * t = (TTree*) f.Get(tree_name); 
  mattak::Waveforms * wf = 0; 
  t->SetBranchAddress(branch_name, &wf); 
  t->GetEntry(i); 
  wf->drawWaveforms(); 
}
