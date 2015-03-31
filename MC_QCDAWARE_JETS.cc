// -*- C++ -*-
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/QCDAware.hh"

#include "Rivet/Analysis.hh"
#include "Rivet/Jet.hh"
#include "Rivet/Projections/FinalPartons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Tools/Logging.hh"

#include "UserInfoParticle.hh"


using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

namespace Rivet {

    class LabeledJet {
        private:
            map<string, Particle> _labelMap;
            Jet _jet;

        public:
            LabeledJet(const Jet& jet) {
                _jet = jet;
                return;
            }

            LabeledJet(const PseudoJet& pj) {
                _jet = Jet(pj);
                return;
            }

            Particle& operator[] (const string& lab) {
                map<string, Particle>::const_iterator p = labelMap.find(lab);

                // automatically make particle if doesn't already
                // exist.
                if (p == labelMap.end())
                    labelMap[lab] = Particle(0, FourMomentum(0, 0, 0, 0));

                return labelMap[lab];
            }

            const Jet jet() {
                return _jet;
            }
    };


    class MC_QCDAWARE_JETS : public Analysis {

        private:
            vector<string> flavors;
            vector<string> labels;
            vector<string> labelsTex;
            vector<string> leadlabs;

            double maxLabelDr;

            QCDAware *qcdawareakt;
            QCDAware *qcdawarekt;

            std::map<string, Histo1DPtr> histos1D;
            std::map<string, Profile1DPtr> profiles1D;
            std::map<string, Histo2DPtr> histos2D;


        public:

            /// Constructor
            // this is really, really ugly.
            MC_QCDAWARE_JETS()
                : Analysis("MC_QCDAWARE_JETS"),
                maxLabelDr(0.2) {

                    flavors.push_back("Unlabeled");
                    flavors.push_back("Gluon");
                    flavors.push_back("Light");
                    flavors.push_back("Charm");
                    flavors.push_back("Bottom");
                    flavors.push_back("Photon");
                    flavors.push_back("Electron");
                    flavors.push_back("Muon");
                    flavors.push_back("Tau");

                    labels.push_back("Akt");
                    labels.push_back("Kt");
                    labels.push_back("MaxPt");
                    labels.push_back("Reclustered");

                    labelsTex.push_back("anti-$k_t$ label");
                    labelsTex.push_back("$k_t$ label");
                    labelsTex.push_back("max-$p_T$ label");
                    labelsTex.push_back("reclustered $k_t$ label");

                    leadlabs.push_back("Inclusive");
                    leadlabs.push_back("Jet0");
                    leadlabs.push_back("Jet1");
                    leadlabs.push_back("Jet2");
                    leadlabs.push_back("Jet3");

                    return;
                }


            void init() {

                FinalPartons fps;
                addProjection(fps, "FinalPartons");

                VisibleFinalState vfs = VisibleFinalState(-2.5, 2.5);
                addProjection(vfs, "VisibleFinalState");

                IdentifiedFinalState lepgammafs(vfs);
                lepgammafs.acceptIdPair(11);
                lepgammafs.acceptIdPair(13);
                lepgammafs.acceptId(22);
                addProjection(lepgammafs, "ElectronsMuonsPhotons");

                TauFinder taufs(TauFinder::ANY);
                addProjection(taufs, "Taus");

                fastjet::contrib::AntiKtMeasure *aktdm =
                    new fastjet::contrib::AntiKtMeasure(0.4);
                fastjet::contrib::KtMeasure *ktdm =
                    new fastjet::contrib::KtMeasure(0.4);

                qcdawareakt = new fastjet::contrib::QCDAware(aktdm);
                qcdawarekt = new fastjet::contrib::QCDAware(ktdm);


                foreach (const string& flav, flavors)
                    foreach (const string& lab, labels)
                    foreach (const string& leadlab, leadlabs)
                    bookLabelHistos(leadlab + "_" + flav + "_" + lab);

                for (unsigned int i = 0; i < labels.size(); i++)
                    for (unsigned int j = i+1; j < labels.size(); j++)
                        foreach (const string& leadlab, leadlabs)
                            bookLabelComparison(leadlab,
                                    labels[i], labelsTex[i],
                                    labels[j], labelsTex[j]);

                return;
            }


            /// Perform the per-event analysis
            void analyze(const Event& event) {
                const double weight = event.weight();

                // first get all final partons
                const Particles& finalPartons =
                    applyProjection<FinalPartons>(event, "FinalPartons").particles();

                const Particles& lepsgammas = 
                    applyProjection<IdentifiedFinalState>(event, "ElectronsMuonsPhotons").particles();

                const Particles& taus =
                    applyProjection<TauFinder>(event, "Taus").taus();

                // first get the qcd-aware parton jets
                vector<Particle> partonJetInputs;

                // loop over final partons
                foreach (const Particle& part, finalPartons) {
                    if (part.abseta() > 7.0)
                        continue;

                    partonJetInputs.push_back(part);
                }

                // leptons and photons
                foreach (const Particle& part, lepsgammas) {
                    // deal with tau decay products later.
                    if (part.fromDecay())
                        continue;

                    if (part.abseta() > 7.0)
                        continue;

                    partonJetInputs.push_back(part);
                }

                // hadronic taus
                foreach (const Particle& tau, taus) {
                    // reject taus from hadron decays and tau chains
                    if (tau.fromDecay())
                        continue;


                    // TODO
                    // this shouldn't be necessary...
                    foreach (const Particle& p, tau.stableDescendants()) {

                        if (tau.abseta() > 7.0)
                            continue;

                        if (p.isHadron()) {
                            partonJetInputs.push_back(tau);
                            break;
                        }
                    }

                    // this isn't a hadronic tau.
                    // all stable descendants of leptonic taus should
                    // be included except neutrinos.
                    foreach (const Particle& part, tau.stableDescendants()) {
                        if (part.isNeutrino())
                            continue;

                        if (part.abseta() > 7.0)
                            continue;

                        partonJetInputs.push_back(part);
                    }
                }


                // make parton jet input pseudojets
                vector<PseudoJet> partonPJs;
                foreach (const Particle& part, partonJetInputs) {
                    PseudoJet partPJ = part.pseudojet();
                    // user_info points to particle.
                    partPJ.set_user_info(new UserInfoParticle(part, "Parton"));

                    // user_index used for flavor-aware clustering.
                    partPJ.set_user_index(part.pid());

                    partonPJs.push_back(partPJ);
                }

                ClusterSequence qcdawareaktcs(partonPJs, qcdawareakt);
                ClusterSequence qcdawarektcs(partonPJs, qcdawarekt);

                const vector<PseudoJet> aktPartonJets =
                    sorted_by_pt(qcdawareaktcs.inclusive_jets(5*GeV));

                const vector<PseudoJet> ktPartonJets =
                    sorted_by_pt(qcdawarektcs.inclusive_jets(5*GeV));

                // now particle jets
                const Particles& visibleParts =
                    applyProjection<VisibleFinalState>(event, "VisibleFinalState").particles();

                // constituents for particle jets
                vector<PseudoJet> particlePJs;
                foreach (const Particle& p, visibleParts) {
                    PseudoJet pj = p.pseudojet();
                    pj.set_user_info(new UserInfoParticle(p, "Particle"));
                    particlePJs.push_back(pj);
                }

                // ghost association of parton jets to particle jets
                foreach (const PseudoJet& aktPJ, aktPartonJets)
                    particlePJs.push_back(
                            ghost(Particle(aktPJ.user_index(), momentum(aktPJ)),
                                "GAAktPartonJet", aktPJ.user_index()));

                foreach (const PseudoJet& ktPJ, ktPartonJets)
                    particlePJs.push_back(
                            ghost(Particle(ktPJ.user_index(), momentum(ktPJ)),
                                "GAKtPartonJet", ktPJ.user_index()));

                // ghost association of final partons to particle jets
                foreach (const Particle& part, partonJetInputs)
                    particlePJs.push_back(ghost(part, "GAFinalParton", part.pid()));

                // ghost association of ALL partons to particle jets
                // for max-pt labeling
                foreach (const GenParticle* gp, Rivet::particles(event.genEvent())) {
                    Particle part(gp);

                    // cut out non-partons and high-eta (including incoming) partons
                    if (!isParton(part) || part.abseta() > 7.0)
                        continue;

                    particlePJs.push_back(ghost(part, "GAParton", part.pid()));
                }

                ClusterSequence akt04cs(particlePJs, JetDefinition(antikt_algorithm, 0.4));

                const vector<PseudoJet> aktJets = sorted_by_pt(akt04cs.inclusive_jets(25*GeV));

                unsigned int iJet = 0;
                foreach (const PseudoJet& j, aktJets) {
                    LabeledJet jet(j);
                    fillJetLabels(jet);

                    fillLabelHistos("Inclusive", jet, weight);

                    switch (iJet) {
                        case 0:
                            fillLabelHistos("Jet0", jet, weight);
                            break;
                        case 1:
                            fillLabelHistos("Jet1", jet, weight);
                            break;
                        case 2:
                            fillLabelHistos("Jet2", jet, weight);
                            break;
                        case 3:
                            fillLabelHistos("Jet3", jet, weight);
                            break;
                    }

                    iJet++;
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


        private:

            void bookLabelHistos(const string& basename) {

                MSG_DEBUG(string("booking label histograms: ") + basename);

                histos1D[basename + "_Pt"] =
                    bookHisto1D(basename + "_Pt", 50, 0, 100*GeV, "$p_T$",
                            "$p_T$ / GeV", "entries");

                histos1D[basename + "_Dpt"] =
                    bookHisto1D(basename + "_Dpt", 50, -1, 1, "$p_T$ resolution",
                            "$p_T$ resolution", "entries");

                histos1D[basename + "_Dr"] =
                    bookHisto1D(basename + "_Dr", 50, 0, 1,
                            "$\\Delta R$", "$\\Delta R$", "entries");

                profiles1D[basename + "_MeanDrVsPt"] =
                    bookProfile1D(basename + "_MeanDrVsPt", 50, 0, 100*GeV,
                            "mean $\\Delta R$ vs $p_T$", "$p_T$ / GeV", "$\\Delta R$");

                profiles1D[basename + "_MeanDptVsDr"] =
                    bookProfile1D(basename + "_MeanDptVsDr", 50, 0, 1,
                            "mean $p_T$ resolution vs $\\Delta R$", "$\\Delta R$", "$p_T$ resolution");

                profiles1D[basename + "_MeanDptVsPt"] =
                    bookProfile1D(basename + "_MeanDptVsPt", 50, 0, 100*GeV,
                            "mean $p_T$ resolution vs $p_T$", "$p_T$ / GeV", "$p_T$ resolution");

                histos2D[basename + "_DrDpt"] = bookHisto2D(basename + "_DrDpt",
                        50, 0, 1, 50, -1, 1,
                        "$\\Delta R$ vs $p_T$ resolution",
                        "\\Delta R", "$p_T$ resolution", "entries");

                return;
            }


            // jet cannot be const because of default Particle return.
            void fillLabelHistos(const string& prefix, LabeledJet& jet, double weight) {

                MSG_DEBUG(string("filling label histograms: ") + prefix);

                double pt = jet.pt();
                cout << "pt: " << pt << endl;
                foreach (const string& lab, labels) {
                    const Particle& labelpart = jet[lab];
                    double dpt = 1 - labelpart.pt()/pt;
                    double dr = deltaR(jet, labelpart);

                    const string basename = prefix + "_" +
                        pidToLabel(labelpart.pid()) + "_" +
                        lab;

                    histos1D[basename + "_Pt"]->fill(pt, weight);
                    histos1D[basename + "_Dpt"]->fill(dpt, weight);
                    histos1D[basename + "_Dr"]->fill(dr, weight);
                    profiles1D[basename + "_MeanDrVsPt"]->fill(pt, dr, weight);
                    profiles1D[basename + "_MeanDptVsDr"]->fill(dr, dpt, weight);
                    profiles1D[basename + "_MeanDptVsPt"]->fill(pt, dpt, weight);
                    histos2D[basename + "_DrDpt"]->fill(dr, dpt, weight);
                }

                for (unsigned int i = 0; i < labels.size(); i++) {
                    for (unsigned int j = i+1; j < labels.size(); j++) {
                        const string basename = prefix + "_" + labels[i] + "LabVs" +
                            labels[j] + "Lab";
                        histos2D[basename]->fill(jet[labels[i]].pid(), jet[labels[j]].pid(), weight);
                    }
                }
            }

            void bookLabelComparison(const string& prefix,
                    const string& lab1, const string& axis1,
                    const string& lab2, const string& axis2) {

                histos2D[prefix + "_" + lab1 + "LabVs" + lab2 + "Lab"] =
                    bookHisto2D(prefix + "_" + lab1 + "LabVs" + lab2 + "Lab",
                            51, -25.5, 25.5, 51, -25.5, 25.5,
                            axis1 + " vs " + axis2,
                            axis1, axis2, "entries");
            }

            // fills in the labels for a given jet
            void fillJetLabels(LabeledJet& jet) {
                MSG_DEBUG("fillng jet labels.");

                FourMomentum jp4 = momentum(jet);

                // for recluster labling
                vector<PseudoJet> partonReclusteredInputs;

                foreach (const PseudoJet& pj, jet.constituents()) {
                    const UserInfoParticle& uip = pj.user_info<UserInfoParticle>();
                    const string& s = uip.str();
                    const Particle& part = uip.particle();

                    if (s == "Particle")
                        continue;


                    // ghost associated partons
                    if (s == "GAFinalParton") {

                        // save pseudojet for reclustering
                        PseudoJet partPJ = part.pseudojet();
                        partPJ.set_user_index(part.pid());
                        partPJ.set_user_info(new UserInfoParticle(part, "Parton"));
                        partonReclusteredInputs.push_back(partPJ);

                        continue;
                    }


                    if (s == "GAParton") {
                        // note the highest-pt parton
                        if (part.pT() > jet["MaxPt"].pT())
                            jet["MaxPt"] = part;

                        continue;
                    }

                    // store best-matched parton label jet
                    if (deltaR(jp4, part) > maxLabelDr)
                        continue;

                    if (s == "GAAktPartonJet" &&
                            (!jet["Akt"].pt() ||
                             deltaR(jp4, part) < deltaR(jp4, jet["Akt"])))
                        jet["Akt"] = part;

                    else if (s == "GAKtPartonJet" &&
                            (!jet["Kt"].pt() ||
                             deltaR(jp4, part) < deltaR(jp4, jet["Kt"])))
                        jet["Kt"] = part;
                }

                // recluster ghost-matched partons
                ClusterSequence qcdawarereclusterktcs(partonReclusteredInputs, qcdawarekt);
                const vector<PseudoJet> reclusterKtPartonJets =
                    sorted_by_pt(qcdawarereclusterktcs.inclusive_jets(5*GeV));

                foreach (const PseudoJet& pj, reclusterKtPartonJets) {
                    if (!jet["Reclustered"].pt() ||
                            deltaR(jp4, momentum(pj)) < deltaR(jp4, jet["Reclustered"]))
                        jet["Reclustered"] = Particle(pj.user_index(), momentum(pj));
                }

                return;
            }


            string pidToLabel(int pid) {

                int abspid = abs(pid);
                switch (abspid) {
                    case 22:
                        return "Photon";
                    case 21:
                        return "Gluon";
                    case 15:
                        return "Tau";
                    case 13:
                        return "Muon";
                    case 11:
                        return "Electron";
                    case 5:
                        return "Bottom";
                    case 4:
                        return "Charm";
                    case 3:
                    case 2:
                    case 1:
                        return "Light";
                    case 0:
                        return "Unlabeled";
                }

                return "";
            }


    };


    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(MC_QCDAWARE_JETS);

}
