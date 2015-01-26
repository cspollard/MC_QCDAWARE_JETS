#include <string>
#include "fastjet/PseudoJet.hh"
#include "Rivet/ParticleBase.hh"

class RivetUserInfo : public fastjet::PseudoJet::UserInfoBase {
    private:
        const Rivet::ParticleBase& _p;
        const std::string& _s;

    public:
        RivetUserInfo(const Rivet::ParticleBase& p, const std::string& s="")
            : _p(p), _s(s) { }

        template <class P>
        const P& get() const {
            return dynamic_cast<const P&>(_p);
        }

        const std::string& str() const {
            return _s;
        }
};


fastjet::PseudoJet ghost(const Rivet::ParticleBase& p, const std::string& s="", const int idx=0) {
    fastjet::PseudoJet pj;
    pj.reset_PtYPhiM(p.pt()*1e-10, p.rap(), p.phi(), p.mass()*1e-10);
    pj.set_user_info(new RivetUserInfo(p, s));
    pj.set_user_index(idx);

    return pj;
}
