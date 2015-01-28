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

                FinalState fsps = FinalState();
                addProjection(fsps, "FinalState");

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




                // book histograms
                hQAktDpt = bookHisto1D("hQAktDpt",
                        50, -1, 1, "Q Akt $p_T$ resolution", "$p_T$ resolution", "entries");
                hQAktDr = bookHisto1D("hQAktDr",
                        50, 0, 1, "Q Akt $\\Delta R / 0.4$", "$\\Delta R / 0.4$", "entries");
                hQKtDpt = bookHisto1D("hQKtDpt",
                        50, -1, 1, "Q Kt $p_T$ resolution", "$p_T$ resolution", "entries");
                hQKtDr = bookHisto1D("hQKtDr",
                        50, 0, 1, "Q Kt $\\Delta R / 0.4$", "$\\Delta R / 0.4$", "entries");
                hQCADpt = bookHisto1D("hQCADpt",
                        50, -1, 1, "Q CA $p_T$ resolution", "$p_T$ resolution", "entries");
                hQCADr = bookHisto1D("hQCADr",
                        50, 0, 1, "Q CA $\\Delta R / 0.4$", "$\\Delta R / 0.4$", "entries");

                hQAktDrDpt = bookHisto2D("hQAktDrDpt",
                        50, 0, 1, 50, -1, 1,
                        "Q Akt $\\Delta R/0.4$ vs $p_T$ resolution", "\\Delta R/0.4", "$p_T$ resolution", "entries");
                hQKtDrDpt = bookHisto2D("hQKtDrDpt",
                        50, 0, 1, 50, -1, 1,
                        "Q Kt $\\Delta R/0.4$ vs $p_T$ resolution", "\\Delta R/0.4", "$p_T$ resolution", "entries");
                hQCADrDpt = bookHisto2D("hQCADrDpt",
                        50, 0, 1, 50, -1, 1,
                        "Q CA $\\Delta R/0.4$ vs $p_T$ resolution", "\\Delta R/0.4", "$p_T$ resolution", "entries");

                hPAktDpt = bookHisto1D("hPAktDpt",
                        50, -1, 1, "P Akt $p_T$ resolution", "$p_T$ resolution", "entries");
                hPAktDr = bookHisto1D("hPAktDr",
                        50, 0, 1, "P Akt $\\Delta R / 0.4$", "$\\Delta R / 0.4$", "entries");
                hPKtDpt = bookHisto1D("hPKtDpt",
                        50, -1, 1, "P Kt $p_T$ resolution", "$p_T$ resolution", "entries");
                hPKtDr = bookHisto1D("hPKtDr",
                        50, 0, 1, "P Kt $\\Delta R / 0.4$", "$\\Delta R / 0.4$", "entries");
                hPCADpt = bookHisto1D("hPCADpt",
                        50, -1, 1, "P CA $p_T$ resolution", "$p_T$ resolution", "entries");
                hPCADr = bookHisto1D("hPCADr",
                        50, 0, 1, "P CA $\\Delta R / 0.4$", "$\\Delta R / 0.4$", "entries");

                hPAktDrDpt = bookHisto2D("hPAktDrDpt",
                        50, 0, 1, 50, -1, 1,
                        "P Akt $\\Delta R/0.4$ vs $p_T$ resolution", "\\Delta R/0.4", "$p_T$ resolution", "entries");
                hPKtDrDpt = bookHisto2D("hPKtDrDpt",
                        50, 0, 1, 50, -1, 1,
                        "P Kt $\\Delta R/0.4$ vs $p_T$ resolution", "\\Delta R/0.4", "$p_T$ resolution", "entries");
                hPCADrDpt = bookHisto2D("hPCADrDpt",
                        50, 0, 1, 50, -1, 1,
                        "P CA $\\Delta R/0.4$ vs $p_T$ resolution", "\\Delta R/0.4", "$p_T$ resolution", "entries");

                hGAktDpt = bookHisto1D("hGAktDpt",
                        50, -1, 1, "G Akt $p_T$ resolution", "$p_T$ resolution", "entries");
                hGAktDr = bookHisto1D("hGAktDr",
                        50, 0, 1, "G Akt $\\Delta R / 0.4$", "$\\Delta R / 0.4$", "entries");
                hGKtDpt = bookHisto1D("hGKtDpt",
                        50, -1, 1, "G Kt $p_T$ resolution", "$p_T$ resolution", "entries");
                hGKtDr = bookHisto1D("hGKtDr",
                        50, 0, 1, "G Kt $\\Delta R / 0.4$", "$\\Delta R / 0.4$", "entries");
                hGCADpt = bookHisto1D("hGCADpt",
                        50, -1, 1, "G CA $p_T$ resolution", "$p_T$ resolution", "entries");
                hGCADr = bookHisto1D("hGCADr",
                        50, 0, 1, "G CA $\\Delta R / 0.4$", "$\\Delta R / 0.4$", "entries");

                hGAktDrDpt = bookHisto2D("hGAktDrDpt",
                        50, 0, 1, 50, -1, 1,
                        "G Akt $\\Delta R/0.4$ vs $p_T$ resolution", "\\Delta R/0.4", "$p_T$ resolution", "entries");
                hGKtDrDpt = bookHisto2D("hGKtDrDpt",
                        50, 0, 1, 50, -1, 1,
                        "G Kt $\\Delta R/0.4$ vs $p_T$ resolution", "\\Delta R/0.4", "$p_T$ resolution", "entries");
                hGCADrDpt = bookHisto2D("hGCADrDpt",
                        50, 0, 1, 50, -1, 1,
                        "G CA $\\Delta R/0.4$ vs $p_T$ resolution", "\\Delta R/0.4", "$p_T$ resolution", "entries");

                hUAktDpt = bookHisto1D("hUAktDpt",
                        50, -1, 1, "U Akt $p_T$ resolution", "$p_T$ resolution", "entries");
                hUAktDr = bookHisto1D("hUAktDr",
                        50, 0, 1, "U Akt $\\Delta R / 0.4$", "$\\Delta R / 0.4$", "entries");
                hUKtDpt = bookHisto1D("hUKtDpt",
                        50, -1, 1, "U Kt $p_T$ resolution", "$p_T$ resolution", "entries");
                hUKtDr = bookHisto1D("hUKtDr",
                        50, 0, 1, "U Kt $\\Delta R / 0.4$", "$\\Delta R / 0.4$", "entries");
                hUCADpt = bookHisto1D("hUCADpt",
                        50, -1, 1, "U CA $p_T$ resolution", "$p_T$ resolution", "entries");
                hUCADr = bookHisto1D("hUCADr",
                        50, 0, 1, "U CA $\\Delta R / 0.4$", "$\\Delta R / 0.4$", "entries");

                hUAktDrDpt = bookHisto2D("hUAktDrDpt",
                        50, 0, 1, 50, -1, 1,
                        "U Akt $\\Delta R/0.4$ vs $p_T$ resolution", "\\Delta R/0.4", "$p_T$ resolution", "entries");
                hUKtDrDpt = bookHisto2D("hUKtDrDpt",
                        50, 0, 1, 50, -1, 1,
                        "U Kt $\\Delta R/0.4$ vs $p_T$ resolution", "\\Delta R/0.4", "$p_T$ resolution", "entries");
                hUCADrDpt = bookHisto2D("hUCADrDpt",
                        50, 0, 1, 50, -1, 1,
                        "U CA $\\Delta R/0.4$ vs $p_T$ resolution", "\\Delta R/0.4", "$p_T$ resolution", "entries");
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

                const Particles& allParts =
                    applyProjection<FinalState>(event, "FinalState").particles();


                // constituents for particle jets
                pjs.clear();
                foreach (const Particle& p, allParts) {
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

                const vector<PseudoJet> aktJets = sorted_by_pt(akt04cs.inclusive_jets(5*GeV));
                const vector<PseudoJet> ktJets = sorted_by_pt(kt04cs.inclusive_jets(5*GeV));
                const vector<PseudoJet> caJets = sorted_by_pt(ca04cs.inclusive_jets(5*GeV));

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
                    }

                    gaKtJets = sortByPt(gaKtJets);
                    gaCAJets = sortByPt(gaCAJets);

                    if (gaAktJets.size()) {
                        gaAktJets = sortByPt(gaAktJets);
                        const Particle& labJet = gaAktJets[0];

                        if (labJet.abspid() <= 6) {
                            hQAktDpt->fill(1 - labJet.pt()/j.pt(), weight);
                            // hQAktDr->fill(deltaR(j.four_mom(), labJet.mom()), weight);
                            // hQAktDrDpt->fill(deltaR(j.four_mom(), labJet), 1 - labJet.pt()/j.pt(), weight);
                        } else if (labJet.abspid() == 22) {
                            hPAktDpt->fill(1 - labJet.pt()/j.pt(), weight);
                            // hPAktDr->fill(deltaR(j.four_mom(), labJet), weight);
                            // hPAktDrDpt->fill(deltaR(j.four_mom(), labJet), 1 - labJet.pt()/j.pt(), weight);
                        } else if (labJet.abspid() == 21) {
                            hGAktDpt->fill(1 - labJet.pt()/j.pt(), weight);
                            // hGAktDr->fill(deltaR(j.four_mom(), labJet), weight);
                            // hGAktDrDpt->fill(deltaR(j.four_mom(), labJet), 1 - labJet.pt()/j.pt(), weight);
                        } else {
                            hUAktDpt->fill(1 - labJet.pt()/j.pt(), weight);
                            // hUAktDr->fill(deltaR(j.four_mom(), labJet), weight);
                            // hUAktDrDpt->fill(deltaR(j.four_mom(), labJet), 1 - labJet.pt()/j.pt(), weight);
                        }
                    }

                    if (gaKtJets.size()) {
                        gaKtJets = sortByPt(gaKtJets);
                        const Particle& labJet = gaKtJets[0];

                        if (labJet.abspid() <= 6) {
                            hQKtDpt->fill(1 - labJet.pt()/j.pt(), weight);
                            // hQKtDr->fill(deltaR(j.four_mom(), labJet.mom()), weight);
                            // hQKtDrDpt->fill(deltaR(j.four_mom(), labJet), 1 - labJet.pt()/j.pt(), weight);
                        } else if (labJet.abspid() == 22) {
                            hPKtDpt->fill(1 - labJet.pt()/j.pt(), weight);
                            // hPKtDr->fill(deltaR(j.four_mom(), labJet), weight);
                            // hPKtDrDpt->fill(deltaR(j.four_mom(), labJet), 1 - labJet.pt()/j.pt(), weight);
                        } else if (labJet.abspid() == 21) {
                            hGKtDpt->fill(1 - labJet.pt()/j.pt(), weight);
                            // hGKtDr->fill(deltaR(j.four_mom(), labJet), weight);
                            // hGKtDrDpt->fill(deltaR(j.four_mom(), labJet), 1 - labJet.pt()/j.pt(), weight);
                        } else {
                            hUKtDpt->fill(1 - labJet.pt()/j.pt(), weight);
                            // hUKtDr->fill(deltaR(j.four_mom(), labJet), weight);
                            // hUKtDrDpt->fill(deltaR(j.four_mom(), labJet), 1 - labJet.pt()/j.pt(), weight);
                        }
                    }

                    if (gaCAJets.size()) {
                        gaCAJets = sortByPt(gaCAJets);
                        const Particle& labJet = gaCAJets[0];

                        if (labJet.abspid() <= 6) {
                            hQCADpt->fill(1 - labJet.pt()/j.pt(), weight);
                            // hQCADr->fill(deltaR(j.four_mom(), labJet.mom()), weight);
                            // hQCADrDpt->fill(deltaR(j.four_mom(), labJet), 1 - labJet.pt()/j.pt(), weight);
                        } else if (labJet.abspid() == 22) {
                            hPCADpt->fill(1 - labJet.pt()/j.pt(), weight);
                            // hPCADr->fill(deltaR(j.four_mom(), labJet), weight);
                            // hPCADrDpt->fill(deltaR(j.four_mom(), labJet), 1 - labJet.pt()/j.pt(), weight);
                        } else if (labJet.abspid() == 21) {
                            hGCADpt->fill(1 - labJet.pt()/j.pt(), weight);
                            // hGCADr->fill(deltaR(j.four_mom(), labJet), weight);
                            // hGCADrDpt->fill(deltaR(j.four_mom(), labJet), 1 - labJet.pt()/j.pt(), weight);
                        } else {
                            hUCADpt->fill(1 - labJet.pt()/j.pt(), weight);
                            // hUCADr->fill(deltaR(j.four_mom(), labJet), weight);
                            // hUCADrDpt->fill(deltaR(j.four_mom(), labJet), 1 - labJet.pt()/j.pt(), weight);
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


            Histo1DPtr hQAktDpt;
            Histo1DPtr hQAktDr;
            Histo1DPtr hQKtDpt;
            Histo1DPtr hQKtDr;
            Histo1DPtr hQCADpt;
            Histo1DPtr hQCADr;

            Histo2DPtr hQAktDrDpt;
            Histo2DPtr hQKtDrDpt;
            Histo2DPtr hQCADrDpt;

            Histo1DPtr hPAktDpt;
            Histo1DPtr hPAktDr;
            Histo1DPtr hPKtDpt;
            Histo1DPtr hPKtDr;
            Histo1DPtr hPCADpt;
            Histo1DPtr hPCADr;

            Histo2DPtr hPAktDrDpt;
            Histo2DPtr hPKtDrDpt;
            Histo2DPtr hPCADrDpt;

            Histo1DPtr hGAktDpt;
            Histo1DPtr hGAktDr;
            Histo1DPtr hGKtDpt;
            Histo1DPtr hGKtDr;
            Histo1DPtr hGCADpt;
            Histo1DPtr hGCADr;

            Histo2DPtr hGAktDrDpt;
            Histo2DPtr hGKtDrDpt;
            Histo2DPtr hGCADrDpt;

            Histo1DPtr hUAktDpt;
            Histo1DPtr hUAktDr;
            Histo1DPtr hUKtDpt;
            Histo1DPtr hUKtDr;
            Histo1DPtr hUCADpt;
            Histo1DPtr hUCADr;

            Histo2DPtr hUAktDrDpt;
            Histo2DPtr hUKtDrDpt;
            Histo2DPtr hUCADrDpt;

    };


    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(MC_QCDAWARE_JETS);

}
