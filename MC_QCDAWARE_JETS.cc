// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalPartons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/contrib/QCDAware.hh"

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

                fastjet::contrib::QCDAware *qcdaware =
                    new fastjet::contrib::QCDAware(new fastjet::contrib::AntiKtMeasure(0.4));

                addProjection(FastJets(fps, qcdaware), "AntiKt04FinalPartonJets");

            }


            /// Perform the per-event analysis
            void analyze(const Event& event) {
                cout << endl << endl << "NEW EVENT" << endl;

                const double weight = event.weight();

                const Particles& partons =
                    applyProjection<FinalPartons>(event, "FinalPartons").particles();

                const Jets& partonJets =
                    applyProjection<FastJets>(event, "AntiKt04FinalPartonJets").jetsByPt();


                cout << "FINAL PARTONS" << endl;
                foreach (const Particle& p, partons) {
                    cout << "pid pt eta : "
                        << p.pid() << " "
                        << p.pt() << " "
                        << p.eta() << endl;
                }

                cout << endl << "FINAL PARTON JETS" << endl;
                foreach (const Jet& j, partonJets) {
                    const fastjet::PseudoJet& p = j.pseudojet();

                    cout << "pid pt eta : "
                        << p.user_index() << " "
                        << p.pt() << " "
                        << p.eta() << endl;
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

            // Data members like post-cuts event weight counters go here


        private:

            /// @name Histograms
            //@{
            Profile1DPtr _h_XXXX;
            Histo1DPtr _h_YYYY;
            //@}


    };



    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(MC_QCDAWARE_JETS);

}
