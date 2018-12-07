set(INCLUDE_LIST
    "${LAB_INCLUDE_ROOT}/FluidsScene.hpp"
#    "${LAB_INCLUDE_ROOT}/Sphere.hpp"
#    "${LAB_INCLUDE_ROOT}/Spring.hpp"
#    "${LAB_INCLUDE_ROOT}/Lagrange.hpp"
    "${LAB_INCLUDE_ROOT}/UniformGrid.hpp"
    "${LAB_INCLUDE_ROOT}/Particle.hpp"
    "${LAB_INCLUDE_ROOT}/ParticleSystem.hpp"
    "${LAB_INCLUDE_ROOT}/MarchingCube.hpp"
    )

set(PATH_INCLUDE "${LAB_INCLUDE_ROOT}/Paths.hpp")
configure_file("${LAB_INCLUDE_ROOT}/Paths.hpp.in" ${PATH_INCLUDE})

set(LAB_INCLUDE_LIST
    ${INCLUDE_LIST}
    ${PATH_INCLUDE}
    PARENT_SCOPE)
