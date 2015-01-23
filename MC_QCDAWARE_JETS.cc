// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/QCDAware.hh"

#include "FinalPartons.hh"


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
            }


            /// Perform the per-event analysis
            void analyze(const Event& event) {
                cout << endl << endl << "NEW EVENT" << endl;

                const double weight = event.weight();

                // first get all final partons
                const Particles& partons =
                    applyProjection<FinalPartons>(event, "FinalPartons").particles();

                // first get the qcd-aware parton jets
                vector<PseudoJet> pjs;
                foreach (const Particle& p, partons) {
                    PseudoJet pj = p.pseudojet();
                    pj.set_user_index(p.pid());
                    pjs.push_back(pj);
                }


                ClusterSequence qcdawareaktcs(pjs, qcdawareakt);
                ClusterSequence qcdawarektcs(pjs, qcdawarekt);
                ClusterSequence qcdawarecacs(pjs, qcdawareca);

                const vector<PseudoJet> aktPartonJets = sorted_by_pt(qcdawareaktcs.inclusive_jets(5*GeV));
                const vector<PseudoJet> ktPartonJets = sorted_by_pt(qcdawarektcs.inclusive_jets(5*GeV));
                const vector<PseudoJet> caPartonJets = sorted_by_pt(qcdawarecacs.inclusive_jets(5*GeV));


                // TODO
                // what's the right way to ghost associate the
                // labeling jets?
                // first add visible particle pseudojets
                // then add parton jet pseudojets w/ negative indices
                // when looking up parton jet indices, we have to
                // subtract nParticles and the size of previous jet
                // collections.
                pjs.clear();

                const Particles& visibleParts =
                    applyProjection<VisibleFinalState>(event, "VisibleFinalState").particles();

                size_t nParts = visibleParts.size();
                for (size_t iPart = 0; iPart < nParts; iPart++) {
                    const Particle& p = visibleParts.at(iPart);
                    PseudoJet pj = p.pseudojet();
                    pj.set_user_index(iPart);
                    pjs.push_back(pj);
                }

                size_t nAktPartJets = aktPartonJets.size();
                for (size_t iJet = 0; iJet < nAktPartJets; iJet++) {
                    PseudoJet pj(aktPartonJets.at(iJet));
                    pj.reset_PtYPhiM(pj.pt()*1e-10, pj.rap(), pj.phi(), pj.m()*1e-10);
                    pj.set_user_index(-iJet);
                    pjs.push_back(pj);
                }

                size_t nKtPartJets = ktPartonJets.size();
                for (size_t iJet = 0; iJet < nKtPartJets; iJet++) {
                    PseudoJet pj(ktPartonJets.at(iJet));
                    pj.set_user_index(-(nAktPartJets+iJet));
                    pj.reset_PtYPhiM(pj.pt()*1e-10, pj.rap(), pj.phi(), pj.m()*1e-10);
                    pjs.push_back(pj);
                }

                size_t nCAPartJets = caPartonJets.size();
                for (size_t iJet = 0; iJet < nCAPartJets; iJet++) {
                    PseudoJet pj(caPartonJets.at(iJet));
                    pj.reset_PtYPhiM(pj.pt()*1e-10, pj.rap(), pj.phi(), pj.m()*1e-10);
                    pj.set_user_index(-(nAktPartJets+nKtPartJets+iJet));
                    pjs.push_back(pj);
                }

                ClusterSequence akt04cs(pjs, JetDefinition(antikt_algorithm, 0.4));
                ClusterSequence kt04cs(pjs, JetDefinition(kt_algorithm, 0.4));
                ClusterSequence ca04cs(pjs, JetDefinition(cambridge_algorithm, 0.4));

                const vector<PseudoJet> aktJets = sorted_by_pt(akt04cs.inclusive_jets(5*GeV));
                const vector<PseudoJet> ktJets = sorted_by_pt(kt04cs.inclusive_jets(5*GeV));
                const vector<PseudoJet> caJets = sorted_by_pt(ca04cs.inclusive_jets(5*GeV));

                cout << endl << "FINAL PARTON AKT JETS" << endl;
                foreach (const PseudoJet& j, aktPartonJets) {
                    cout << "pid pt eta phi: "
                        << j.user_index() << " "
                        << j.pt() << " "
                        << j.eta() << " "
                        << j.phi() << endl;
                }

                cout << endl << "FINAL PARTICLE AKT JETS" << endl;
                foreach (const PseudoJet& j, aktJets) {
                    cout << "pt eta phi: "
                        << j.pt() << " "
                        << j.eta() << " "
                        << j.phi() << endl;

                    cout << "associated parton jets:" << endl;
                    foreach (const PseudoJet& pj, j.constituents()) {
                        // not a ghost-associated pseudojet
                        if (pj.user_index() >= 0) continue;

                        int idx = abs(pj.user_index());

                        // TODO
                        // this is stupid.
                        PseudoJet labjet;
                        if (idx >= nAktPartJets + nKtPartJets) {
                            labjet = caPartonJets.at(idx - nAktPartJets - nKtPartJets);
                            cout << "\tcamkt pid pt eta phi: "
                                << labjet.user_index() << " "
                                << labjet.pt() << " "
                                << labjet.eta() << " "
                                << labjet.phi() << endl;
                        } else if (idx >= nAktPartJets) {
                            labjet = ktPartonJets.at(idx - nAktPartJets);
                            cout << "\tkt pid pt eta phi: "
                                << labjet.user_index() << " "
                                << labjet.pt() << " "
                                << labjet.eta() << " "
                                << labjet.phi() << endl;
                        } else if (idx >= 0) {
                            labjet = aktPartonJets.at(idx);
                            cout << "\takt pid pt eta phi: "
                                << labjet.user_index() << " "
                                << labjet.pt() << " "
                                << labjet.eta() << " "
                                << labjet.phi() << endl;
                        } else {
                            cout << "ERROR!!! This should be a ghost particle, but it isnt?!" << endl;
                        }
                    }
                }


                return;
            }


            /// Normalise histograms etc., after the run
            void finalize() {

                /// @todo Normalise, scale and otherwise manipulate histograms here

                // scale(_h_YYYY, crossSection()/sumOfWeights()); // norm to cross section
                // normalize(_h_YYYY); // normalize to unity

            }

            //@}


        private:
            QCDAware *qcdawareakt;
            QCDAware *qcdawarekt;
            QCDAware *qcdawareca;


    };



    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(MC_QCDAWARE_JETS);

}
