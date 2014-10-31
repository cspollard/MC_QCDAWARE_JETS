// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalPartons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/QCDAware.hh"


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

                fastjet::contrib::QCDAwareDistanceMeasure<AntiKtMeasure> *dm =
                    new fastjet::contrib::QCDAwareDistanceMeasure<AntiKtMeasure>(0.4);
                qcdaware = new fastjet::contrib::QCDAware(dm);
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
                    // TODO
                    // TEST
                    if (abs(id) < 7 || id == 21 || id == 22)
                        id = 0;
                    // if (id == 21 || id == 22)
                        // id = 0;

                    pj.set_user_index(id);
                    pjs.push_back(pj);
                }


                ClusterSequence qcdawarecs(pjs, qcdaware);
                const vector<PseudoJet> partonJets = sorted_by_pt(qcdawarecs.inclusive_jets(25*GeV));

                cout << "FINAL PARTONS" << endl;
                foreach (const Particle& p, partons) {
                    cout << "pid pt eta phi : "
                        << p.pid() << " "
                        << p.pt() << " "
                        << p.eta() << " "
                        << p.phi() << endl;
                }

                cout << endl << "FINAL PARTON JETS" << endl;
                foreach (const PseudoJet& j, partonJets) {
                    cout << "pid pt eta phi : "
                        << j.user_index() << " "
                        << j.pt() << " "
                        << j.eta() << " "
                        << j.phi() << endl;
                }

                const Jets& testjets =
                    applyProjection<FastJets>(event, "AntiKt04FinalPartonJets").jetsByPt(25*GeV);
                cout << endl << "FINAL PARTON TEST JETS" << endl;
                foreach (const Jet& j, testjets) {
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
            QCDAware *qcdaware;


    };



    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(MC_QCDAWARE_JETS);

}
