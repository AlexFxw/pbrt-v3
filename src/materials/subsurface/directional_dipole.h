#ifndef PBRT_V3_DIRECTIONAL_DIPOLE_H
#define PBRT_V3_DIRECTIONAL_DIPOLE_H

#include "pbrt.h"
#include "bssrdf.h"

namespace pbrt {

class DirectionalBSSRDF : public SeparableBSSRDF {
public:
    bool ContainSingleScattering() const override { return true; }

    Spectrum Sr(Float d) const;
    Float Sample_Sr(int ch, Float u) const;
    Float Pdf_Sr(int ch, Float r) const;

private:
};

}

#endif //PBRT_V3_DIRECTIONAL_DIPOLE_H
