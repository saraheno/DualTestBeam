#ifndef DUALCRYS_CALORIMETER_CONVERTER_H
#define DUALCRYS_CALORIMETER_CONVERTER_H

#include "Gaudi/Property.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "DualCrysCalorimeterHit.h"
#include "DualCrysCalDigi.h"

class DualCrysCalConverter {
public:
    DualCrysCalConverter(float threshold, float calibrationCoeff)
        : m_threshold(threshold), m_calibrationCoeff(calibrationCoeff) {}

    edm4hep::CalorimeterHitCollection convert(const std::vector<CalVision::DualCrysCalorimeterHit>& simHits) {
        edm4hep::CalorimeterHitCollection recoHits;
        
        for (const auto& hit : simHits) {
            float energy = hit.edeprelativistic * m_calibrationCoeff;
            if (energy < m_threshold) continue; // Apply threshold filter
            
            edm4hep::CalorimeterHit recoHit;
            recoHit.setEnergy(energy);
            recoHit.setPosition({hit.xmax, hit.ymax, 0});
            recoHits.push_back(recoHit);
        }
        
        return recoHits;
    }

private:
    float m_threshold;
    float m_calibrationCoeff;
};

#endif // DUALCRYS_CALORIMETER_CONVERTER_H
