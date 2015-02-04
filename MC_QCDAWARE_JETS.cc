// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/QCDAware.hh"

#include "FinalPartons.hh"

#include "UserInfoParticle.hh"


using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

namespace Rivet {

    class MC_QCDAWARE_JETS : public Analysis {
        public:

            /// Constructor
            MC_QCDAWARE_JETS()
                : Analysis("MC_QCDAWARE_JETS")
            {    }


        public:

            /// @name Analysis methods
            //@{

            void init() {

                FinalPartons fps = FinalPartons();
                addProjection(fps, "FinalPartons");

                VisibleFinalState vfs = VisibleFinalState();
                addProjection(vfs, "VisibleFinalState");

                addProjection(FastJets(fps, FastJets::ANTIKT, 0.4), "AntiKt04FinalPartonJets");
                addProjection(FastJets(fps, FastJets::KT, 0.4), "Kt04FinalPartonJets");
                addProjection(FastJets(fps, FastJets::CAM, 0.4), "CA04FinalPartonJets");

                fastjet::contrib::QCDAwareDistanceMeasure<AntiKtMeasure> *aktdm =
                    new fastjet::contrib::QCDAwareDistanceMeasure<AntiKtMeasure>(0.4);
                fastjet::contrib::QCDAwareDistanceMeasure<KtMeasure> *ktdm =
                    new fastjet::contrib::QCDAwareDistanceMeasure<KtMeasure>(0.4);
                fastjet::contrib::QCDAwareDistanceMeasure<CAMeasure> *cadm =
                    new fastjet::contrib::QCDAwareDistanceMeasure<CAMeasure>(0.4);

                qcdawareakt = new fastjet::contrib::QCDAware(aktdm);
                qcdawarekt = new fastjet::contrib::QCDAware(ktdm);
                qcdawareca = new fastjet::contrib::QCDAware(cadm);




                // book labeling histograms
                bookLabelHistos("GluonAkt");
                bookLabelHistos("LightAkt");
                bookLabelHistos("CharmAkt");
                bookLabelHistos("BottomAkt");
                bookLabelHistos("PhotonAkt");
                bookLabelHistos("MuonAkt");
                bookLabelHistos("ElectronAkt");

                bookLabelHistos("GluonKt");
                bookLabelHistos("LightKt");
                bookLabelHistos("CharmKt");
                bookLabelHistos("BottomKt");
                bookLabelHistos("PhotonKt");
                bookLabelHistos("MuonKt");
                bookLabelHistos("ElectronKt");

                // book unlabeled Pt histogram
                histos1D["UnlabeledAktPt"] = bookHisto1D("UnlabeledAktPt",
                        50, 0, 100*GeV, "$p_T$", "$p_T$ / GeV", "entries");
                histos1D["UnlabeledKtPt"] = bookHisto1D("UnlabeledKtPt",
                        50, 0, 100*GeV, "$p_T$", "$p_T$ / GeV", "entries");
            }


            /// Perform the per-event analysis
            void analyze(const Event& event) {
                const double weight = event.weight();

                // first get all final partons
                const Particles& partons =
                    applyProjection<FinalPartons>(event, "FinalPartons").particles();

                // first get the qcd-aware parton jets
                vector<PseudoJet> pjs;
                PseudoJet tmpPJ;
                foreach (const Particle& parton, partons) {
                    tmpPJ = parton.pseudojet();
                    // user_info points to particle.
                    tmpPJ.set_user_info(new UserInfoParticle(parton));

                    // user_index used for flavor-aware clustering.
                    tmpPJ.set_user_index(parton.pid());

                    pjs.push_back(tmpPJ);
                }


                ClusterSequence qcdawareaktcs(pjs, qcdawareakt);
                ClusterSequence qcdawarektcs(pjs, qcdawarekt);
                ClusterSequence qcdawarecacs(pjs, qcdawareca);

                const vector<PseudoJet> aktPartonJets =
                    sorted_by_pt(qcdawareaktcs.inclusive_jets(5*GeV));

                const vector<PseudoJet> ktPartonJets =
                    sorted_by_pt(qcdawarektcs.inclusive_jets(5*GeV));

                const vector<PseudoJet> caPartonJets =
                    sorted_by_pt(qcdawarecacs.inclusive_jets(5*GeV));


                // now particle jets
                const Particles& visibleParts =
                    applyProjection<VisibleFinalState>(event, "VisibleFinalState").particles();

                // constituents for particle jets
                pjs.clear();
                foreach (const Particle& p, visibleParts) {
                    tmpPJ = p.pseudojet();
                    tmpPJ.set_user_info(new UserInfoParticle(p));
                    pjs.push_back(tmpPJ);
                }

                // ghost association of parton jets to particle jets
                foreach (const PseudoJet& aktPJ, aktPartonJets)
                    pjs.push_back(ghost(Particle(aktPJ.user_index(), momentum(aktPJ)),
                                "GhostAktPartonJet"));

                foreach (const PseudoJet& ktPJ, ktPartonJets)
                    pjs.push_back(ghost(Particle(ktPJ.user_index(), momentum(ktPJ)),
                                "GhostKtPartonJet"));

                foreach (const PseudoJet& caPJ, caPartonJets)
                    pjs.push_back(ghost(Particle(caPJ.user_index(), momentum(caPJ)),
                                "GhostCAPartonJet"));

                ClusterSequence akt04cs(pjs, JetDefinition(antikt_algorithm, 0.4));
                ClusterSequence kt04cs(pjs, JetDefinition(kt_algorithm, 0.4));
                ClusterSequence ca04cs(pjs, JetDefinition(cambridge_algorithm, 0.4));

                const vector<PseudoJet> aktJets = sorted_by_pt(akt04cs.inclusive_jets(25*GeV));
                const vector<PseudoJet> ktJets = sorted_by_pt(kt04cs.inclusive_jets(25*GeV));
                const vector<PseudoJet> caJets = sorted_by_pt(ca04cs.inclusive_jets(25*GeV));

                foreach (const PseudoJet& j, aktJets) {

                    Particles gaAktJets;
                    Particles gaKtJets;
                    Particles gaCAJets;

                    foreach (const PseudoJet& p, j.constituents()) {
                        const UserInfoParticle& uip = p.user_info<UserInfoParticle>();
                        const string& s = uip.str();

                        // not a ghost-associated pseudojet
                        if (s == "") continue;

                        const Particle& tmpPart = uip.particle();

                        if (s == "GhostAktPartonJet")
                            gaAktJets.push_back(tmpPart);

                        else if (s == "GhostKtPartonJet")
                            gaKtJets.push_back(tmpPart);

                        else if (s == "GhostCAPartonJet")
                            gaCAJets.push_back(tmpPart);

                        else
                            cout << "ERROR: unrecognized jet label!!!!!" << endl;
                    }

                    gaAktJets = sortByPt(gaAktJets);
                    gaKtJets = sortByPt(gaKtJets);
                    gaCAJets = sortByPt(gaCAJets);

                    if (gaAktJets.size()) {
                        const Particle& labJet = gaAktJets[0];

                        string name = pidToLabel(labJet.pid()) + "Akt";
                        fillLabelHistos(name, weight, Jet(j).mom(), labJet.mom());
                    } else
                        histos1D["UnlabeledAktPt"]->fill(j.pt(), weight);

                    if (gaKtJets.size()) {
                        const Particle& labJet = gaKtJets[0];

                        string name = pidToLabel(labJet.pid()) + "Kt";
                        fillLabelHistos(name, weight, Jet(j).mom(), labJet.mom());
                    } else
                        histos1D["UnlabeledKtPt"]->fill(j.pt(), weight);

                }

                return;
            }


            /// Normalise histograms etc., after the run
            void finalize() {

                // normalize to 1/fb
                double norm = 1000*crossSection()/sumOfWeights();
                for (map< string, Histo1DPtr>::iterator p = histos1D.begin(); p != histos1D.end(); ++p)
                    p->second->scaleW(norm); // norm to cross section

                for (map< string, Profile1DPtr>::iterator p = profiles1D.begin(); p != profiles1D.end(); ++p)
                    p->second->scaleW(norm); // norm to cross section

                for (map< string, Histo2DPtr>::iterator p = histos2D.begin(); p != histos2D.end(); ++p)
                    p->second->scaleW(norm); // norm to cross section

            }

            //@}


        private:
            QCDAware *qcdawareakt;
            QCDAware *qcdawarekt;
            QCDAware *qcdawareca;

            std::map<string, Histo1DPtr> histos1D;
            std::map<string, Profile1DPtr> profiles1D;
            std::map<string, Histo2DPtr> histos2D;


            void bookLabelHistos(const string& basename) {

                histos1D[basename + "Pt"] =
                    bookHisto1D(basename + "Pt", 50, 0, 100*GeV, "$p_T$",
                            "$p_T$ / GeV", "entries");

                histos1D[basename + "Dpt"] =
                    bookHisto1D(basename + "Dpt", 50, -1, 1, "$p_T$ resolution",
                            "$p_T$ resolution", "entries");

                histos1D[basename + "Dr"] =
                    bookHisto1D(basename + "Dr", 50, 0, 1,
                            "$\\Delta R$", "$\\Delta R$", "entries");

                profiles1D[basename + "MeanDrVsPt"] =
                    bookProfile1D(basename + "MeanDrVsPt", 50, 0, 100*GeV,
                            "mean $\\Delta R$ vs $p_T$", "$p_T$ / GeV", "$\\Delta R$");

                profiles1D[basename + "MeanDptVsDr"] =
                    bookProfile1D(basename + "MeanDptVsDr", 50, 0, 1,
                            "mean $p_T$ resolution vs $\\Delta R$", "$\\Delta R$", "$p_T$ resolution");

                profiles1D[basename + "MeanDptVsPt"] =
                    bookProfile1D(basename + "MeanDptVsPt", 50, 0, 100*GeV,
                            "mean $p_T$ resolution vs $p_T$", "$p_T$ / GeV", "$p_T$ resolution");

                histos2D[basename + "DrDpt"] = bookHisto2D(basename + "DrDpt",
                        50, 0, 1, 50, -1, 1,
                        "$\\Delta R$ vs $p_T$ resolution",
                        "\\Delta R", "$p_T$ resolution", "entries");

                return;
            }


            void fillLabelHistos(const string& basename, double weight,
                    const FourMomentum& jet, const FourMomentum& label) {

                double pt = jet.pt();
                double dpt = 1 - label.pt()/pt;
                double dr = deltaR(jet, label);

                histos1D[basename + "Pt"]->fill(pt, weight);
                histos1D[basename + "Dpt"]->fill(dpt, weight);
                histos1D[basename + "Dr"]->fill(dr, weight);
                profiles1D[basename + "MeanDrVsPt"]->fill(pt, dr, weight);
                profiles1D[basename + "MeanDptVsDr"]->fill(dr, dpt, weight);
                profiles1D[basename + "MeanDptVsPt"]->fill(pt, dpt, weight);
                histos2D[basename + "DrDpt"]->fill(dr, dpt, weight);

            }


            string pidToLabel(int pid) {

                int abspid = abs(pid);
                string name;
                switch (abspid) {
                    case 22:
                        name = "Photon";
                        break;
                    case 21:
                        name = "Gluon";
                        break;
                    case 13:
                        name = "Muon";
                        break;
                    case 11:
                        name = "Electron";
                        break;
                    case 5:
                        name = "Bottom";
                        break;
                    case 4:
                        name = "Charm";
                        break;
                    case 3:
                    case 2:
                    case 1:
                        name = "Light";
                        break;
                }

                return name;
            }


    };


    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(MC_QCDAWARE_JETS);

}
