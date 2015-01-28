// -*- C++ -*-

#include "FinalPartons.hh"
#include "Rivet/Tools/RivetHepMC.hh"

namespace Rivet {

    bool FinalPartons::accept(const Particle& p) const {

        if (!isParton(p) && !isPhoton(p))
            return false;

        // only getting gluons for some reason...
        foreach (const Particle& c, p.children())
            if (isParton(c) || isPhoton(c))
                return false;

        return _cuts->accept(p);
    }


    void FinalPartons::project(const Event& e) {
        _theParticles.clear();

        foreach (const GenParticle* gp, Rivet::particles(e.genEvent())) {
            if (!gp)
                continue;

            const Particle& p = Particle(gp);

            if (accept(p))
                _theParticles.push_back(p);
        }

        return;
    }

}
