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

                const Particles& partons =
                    applyProjection<FinalPartons>(event, "FinalPartons").particles();

                vector<PseudoJet> pjs;
                foreach (const Particle& p, partons) {
                    PseudoJet pj = p;
                    int id = p.pid();
                    if (id == 21 || id == 22)
                        id = 0;

                    pj.set_user_index(id);
                    pjs.push_back(pj);
                }


                ClusterSequence qcdawareaktcs(pjs, qcdawareakt);
                ClusterSequence qcdawarektcs(pjs, qcdawarekt);
                ClusterSequence qcdawarecacs(pjs, qcdawareca);

                const vector<PseudoJet> aktPartonJets = sorted_by_pt(qcdawareaktcs.inclusive_jets(5*GeV));
                const vector<PseudoJet> ktPartonJets = sorted_by_pt(qcdawarektcs.inclusive_jets(5*GeV));
                const vector<PseudoJet> caPartonJets = sorted_by_pt(qcdawarecacs.inclusive_jets(5*GeV));

                cout << "FINAL PARTONS" << endl;
                foreach (const Particle& p, sortByPt(partons)) {
                    cout << "pid pt eta phi : "
                        << p.pid() << " "
                        << p.pt() << " "
                        << p.eta() << " "
                        << p.phi() << endl;
                }

                cout << endl << "FINAL PARTON AKT JETS" << endl;
                foreach (const PseudoJet& j, aktPartonJets) {
                    cout << "pid pt eta phi : "
                        << j.user_index() << " "
                        << j.pt() << " "
                        << j.eta() << " "
                        << j.phi() << endl;
                }

                cout << endl << "FINAL PARTON AKT TEST JETS" << endl;
                foreach (const Jet& j,
                        applyProjection<FastJets>(event, "AntiKt04FinalPartonJets").jetsByPt(5*GeV)) {
                    cout << "pt eta phi : "
                        << j.pt() << " "
                        << j.eta() << " "
                        << j.phi() << endl;
                }

                cout << endl << "FINAL PARTON KT JETS" << endl;
                foreach (const PseudoJet& j, ktPartonJets) {
                    cout << "pid pt eta phi : "
                        << j.user_index() << " "
                        << j.pt() << " "
                        << j.eta() << " "
                        << j.phi() << endl;
                }

                cout << endl << "FINAL PARTON KT TEST JETS" << endl;
                foreach (const Jet& j,
                        applyProjection<FastJets>(event, "Kt04FinalPartonJets").jetsByPt(5*GeV)) {
                    cout << "pt eta phi : "
                        << j.pt() << " "
                        << j.eta() << " "
                        << j.phi() << endl;
                }

                cout << endl << "FINAL PARTON CA JETS" << endl;
                foreach (const PseudoJet& j, caPartonJets) {
                    cout << "pid pt eta phi : "
                        << j.user_index() << " "
                        << j.pt() << " "
                        << j.eta() << " "
                        << j.phi() << endl;
                }

                cout << endl << "FINAL PARTON CA TEST JETS" << endl;
                foreach (const Jet& j,
                        applyProjection<FastJets>(event, "CA04FinalPartonJets").jetsByPt(5*GeV)) {
                    cout << "pt eta phi : "
                        << j.pt() << " "
                        << j.eta() << " "
                        << j.phi() << endl;
                }


                /// @todo Do the event by event analysis here

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
