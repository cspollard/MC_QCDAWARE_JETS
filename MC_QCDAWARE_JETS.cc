// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/QCDAware.hh"

#include "FinalPartons.hh"

#include "RivetUserInfo.hh"


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
                    // user_info points to particle.
                    pj.set_user_info(new RivetUserInfo(p));

                    // user_index used for flavor-aware clustering.
                    pj.set_user_index(p.pid());

                    pjs.push_back(pj);
                }


                ClusterSequence qcdawareaktcs(pjs, qcdawareakt);
                ClusterSequence qcdawarektcs(pjs, qcdawarekt);
                ClusterSequence qcdawarecacs(pjs, qcdawareca);

                const vector<PseudoJet> aktPartonJets = sorted_by_pt(qcdawareaktcs.inclusive_jets(5*GeV));
                const vector<PseudoJet> ktPartonJets = sorted_by_pt(qcdawarektcs.inclusive_jets(5*GeV));
                const vector<PseudoJet> caPartonJets = sorted_by_pt(qcdawarecacs.inclusive_jets(5*GeV));


                const Particles& visibleParts =
                    applyProjection<VisibleFinalState>(event, "VisibleFinalState").particles();


                // constituents for particle jets
                pjs.clear();
                foreach (const Particle& p, visibleParts) {
                    PseudoJet pj = p.pseudojet();
                    pj.set_user_info(new RivetUserInfo(p));
                    pjs.push_back(pj);
                }

                // ghost association
                foreach (const PseudoJet& pj, aktPartonJets) {
                    const Particle p(pj.user_index(), momentum(pj));
                    pjs.push_back(ghost(p, "GhostAktPartonJet"));
                }

                foreach (const PseudoJet& pj, ktPartonJets) {
                    const Particle p(pj.user_index(), momentum(pj));
                    pjs.push_back(ghost(p, "GhostKtPartonJet"));
                }

                foreach (const PseudoJet& pj, ktPartonJets) {
                    const Particle p(pj.user_index(), momentum(pj));
                    pjs.push_back(ghost(p, "GhostCAPartonJet"));
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

                cout << endl << "FINAL PARTON KT JETS" << endl;
                foreach (const PseudoJet& j, ktPartonJets) {
                    cout << "pid pt eta phi: "
                        << j.user_index() << " "
                        << j.pt() << " "
                        << j.eta() << " "
                        << j.phi() << endl;
                }

                cout << endl << "FINAL PARTON CA JETS" << endl;
                foreach (const PseudoJet& j, caPartonJets) {
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
                        const RivetUserInfo& rui =
                            pj.user_info<RivetUserInfo>();

                        // not a ghost-associated pseudojet
                        if (rui.str() == "") continue;

                        const Particle& labjet = rui.get<Particle>();
                        cout << "\t" << rui.str() << endl
                            << "\tpid pt eta phi: "
                            << labjet.pid() << " "
                            << labjet.pt() << " "
                            << labjet.eta() << " "
                            << labjet.phi() << endl;
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
