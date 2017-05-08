#include "material.h"

Material::Material(MaterialType mt) : materialType(mt)
{
}

std::string Material::getMaterialName()
{
    switch(materialType) {

    case CONCRETE:
        return "Concrete";
    case STRUCTURAL_STEEL:
        return "Structural Steel";
    case REINFORCEMENT:
        return "Reinforcement";
    case REINFORCEMENT_BARS:
        return "Reinforcement Bars";
    case FAFITIS_CONCRETE:
        return "Fafitis Concrete";
    case FAFITIS_REINFORCEMENT:
        return "Fafitis Reinforcement";
    case HOLLOW:
        return "Hollow";
    case FIBER_REINFORCED_PLASTIC:
        return "FPR";
    default:
        return "No Data";
    }
}

MaterialType Material::getMaterialType() {
    return materialType;
}
