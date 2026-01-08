
/// simple program to get a quick estimate of cross sections and a
/// handful of basic distributions from MC run and flat/HepMC/UW
/// files, etc.
///
#include "AnalysisFramework.hh"
#include "boost/foreach.hpp"
#define foreach BOOST_FOREACH


using namespace std;
using namespace fastjet;

double log10(double x) {return log(x)/log(10.0);}

/// Example class derived from AnalysisFramework that will help in evaluation of
/// some basic cross sections. 
class XSctAnalysis : public AnalysisFramework {
public:

  Selector jet_selector;
  bool alexander_orig_analysis;
  
  XSctAnalysis(CmdLine & cmdline) : AnalysisFramework(cmdline) {}
  void user_startup() {
    // extra user parameters
    param["jj.deltay.min"] = 4.5;
    param["jj.m.min"] = 600.0;
    param["jet.ptmin"] = 25.0;
    param["jet.rapmax"] = 4.5;

    jet_def = JetDefinition(antikt_algorithm, cmdline.value<double>("-R",0.4));

    DefaultHist::set_defaults(0.0, 4.4, cmdline.value("-hist.binwidth",0.2));

    alexander_orig_analysis = cmdline.present("-alexander");

    if (alexander_orig_analysis) {
      jet_selector = SelectorPtMin(param["jet.ptmin"]);
    } else {
      jet_selector = SelectorPtMin(param["jet.ptmin"]) && SelectorAbsRapMax(param["jet.rapmax"]);
    }
  }

  void user_post_startup() {
    header << "# jet_selector = " << jet_selector.description() << endl;
  }
  
  void analyse_event() {
    // cout << "Cross section pointer: " 
    // 	 << driver->generator->gen_event()->cross_section()  << " "
    // 	 << driver->generator->gen_event()->weights().size()
    // 	 << endl;

    double evwgt = driver->generator->hadron_level().weight();
    xsections["total cross section"] += evwgt;
    gen_hists["total.xsct"].set_lims_add_entry(0.0, 1.0, 1.0, 0.5, evwgt);
    vector<PseudoJet> particles = driver->generator->hadron_level().particles();

    // takes things from a Higgs or a Z (in case we do VBF -> Z)
    Selector sel_from_higgs = SelectorPDGId(25) || SelectorPDGId(23) ||
        SelectorComesFrom(driver->generator->gen_event(), SelectorPDGId(25) || SelectorPDGId(23));
    vector<PseudoJet> hadronic, higgs_vect;
    sel_from_higgs.sift(particles, higgs_vect, hadronic);
    PseudoJet higgses = join(higgs_vect);

    norm_hists["higgs.mHH"].set_lims_add_entry(0.0, 200.0, 2.0, higgses.m(), evwgt);
    gen_hists["incl.ptHH"].set_lims_add_entry(0.0, 600.0, 5.0, higgses.pt(), evwgt);
    

    // VBF cuts
    vector<PseudoJet> jets = jet_selector(jet_def(hadronic));
    bool pass_cuts =
                (jets.size() >= 2U) &&
                abs(jets[0].rap() - jets[1].rap()) >= param["jj.deltay.min"] &&
                (jets[0]+jets[1]).m() > param["jj.m.min"] &&
                (jets[0].rap() * jets[1].rap()) < 0;
    if (pass_cuts && alexander_orig_analysis) {
      pass_cuts =  abs(jets[0].rap()) < param["jet.rapmax"]
                && abs(jets[1].rap()) < param["jet.rapmax"];
    }
    gen_hists["vbfcut.xsct"].set_lims_add_entry(0.0, 1.0, 1.0, pass_cuts? 0.5 : -0.5, evwgt);
    if (!pass_cuts) return;
    xsections["VBF-cut cross section"] += evwgt;

    if (jets.size() >= 3U) xsections["VBF-cut and >= 3jets"] += evwgt;
    if (jets.size() >= 4U) xsections["VBF-cut and >= 4jets"] += evwgt;
    
    // now get some histograms
    gen_hists["vbfcut.ptj1"].set_lims_add_entry(0.0, 600.0,  5.0, jets[0].pt(), evwgt);    
    gen_hists["vbfcut.ptj2"].set_lims_add_entry(0.0, 600.0,  5.0, jets[1].pt(), evwgt);
    gen_hists["vbfcut.ptHH" ].set_lims_add_entry(0.0, 600.0,  5.0,   higgses.pt(), evwgt);
    
    gen_hists["vbfcut.ptj2"].set_lims_add_entry(0.0, 600.0,  5.0, jets[1].pt(), evwgt);
    
    
    gen_hists["vbfcut.yj1" ].set_lims_add_entry(-6.0,  6.0,  0.1, jets[0].rap(), evwgt);
    gen_hists["vbfcut.yj2" ].set_lims_add_entry(-6.0,  6.0,  0.1, jets[1].rap(), evwgt);    
    gen_hists["vbfcut.yHH"  ].set_lims_add_entry(-6.0,  6.0,  0.1,   higgses.rap(), evwgt);
    
    gen_hists["vbfcut.dyj1j2"].set_lims_add_entry(-10.0, 10.0, 0.5, jets[1].rap()-jets[0].rap(), evwgt);
    gen_hists["vbfcut.absdyj1j2"].set_lims_add_entry(0, 10.0, 0.5, abs(jets[1].rap()-jets[0].rap()), evwgt);
    gen_hists["vbfcut.dphij1j2"].set_lims_add_entry(0, pi, pi/12, abs(jets[1].delta_phi_to(jets[0])), evwgt);
    gen_hists["vbfcut.dRj1j2"].set_lims_add_entry(4.0, 10.0, 0.5, jets[1].delta_R(jets[0]), evwgt);
  
    gen_hists["vbfcut.Mjj" ].set_lims_add_entry(0.0, 6000.0, 50.0, (jets[0]+jets[1]).m(), evwgt);     

    double htjets = SelectorIdentity().scalar_pt_sum(jets);
    gen_hists["vbfcut.HTjets"].set_lims_add_entry(0.0, 1000.0, 10.0, htjets, evwgt);
    gen_hists["vbfcut.HTall"].set_lims_add_entry(0.0, 1000.0, 10.0, htjets + higgses.pt(), evwgt);
    
    if (jets.size() >=3U) {
      gen_hists["vbfcut.ptj3"].set_lims_add_entry(0.0, 600.0,  5.0, jets[2].pt(), evwgt);
      gen_hists["vbfcut.yj3" ].set_lims_add_entry(-6.0,  6.0,  0.1, jets[2].rap(), evwgt);
      const PseudoJet * refjet;
      if (abs(jets[0].rap()-jets[2].rap()) < abs(jets[1].rap()-jets[2].rap())) {
        refjet = & jets[0];
      } else {
        refjet = & jets[1];
      }
      // define this variable such that positive means central and negative means forward of tagging jets
      double minrap_j1j3_j2j3 = jets[2].rap() - refjet->rap();
      if (refjet->rap() > 0) minrap_j1j3_j2j3 *= -1;

      //cout << jets[0].rap() << " " << jets[1].rap() << " " << refjet->rap() << " " << jets[2].rap() << " " << minrap_j1j3_j2j3 << endl;
      
      double minR_j1j3_j2j3 = refjet->delta_R(jets[2]);
      if (minrap_j1j3_j2j3 < 0) minR_j1j3_j2j3 *= -1;
      
      gen_hists["vbfcut.minrap_j1j2_j2j3" ].set_lims_add_entry(-6.0,  6.0,  0.1, minrap_j1j3_j2j3, evwgt);
      gen_hists["vbfcut.minR_j1j2_j2j3" ].set_lims_add_entry(-6.0,  6.0,  0.1, minR_j1j3_j2j3, evwgt);
      gen_hists["vbfcut.ystarj3"].set_lims_add_entry(-6.0, 6.0, 0.1, jets[2].rap() - 0.5*(jets[0].rap()+jets[1].rap()), evwgt);
    }
    
  }
};

//----------------------------------------------------------------------
int main (int argc, char ** argv) {
  
  CmdLine cmdline(argc,argv);
  XSctAnalysis analysis(cmdline);
  analysis.run();
}
