#pragma once

#include <atlas/utils/Geometry.hpp>
#include <atlas/gl/Buffer.hpp>
#include <atlas/gl/VertexArrayObject.hpp>
#include <cmath>
#include <map>

#include "Particle.hpp"
#include "UniformGrid.hpp"


#define PI 3.1415926535897932384626422832795028841971
//#define EPSILON 0.1
#define EPSILON 0.0001

#define ITERATIONS 4
#define PARTICLE_DIM 8
namespace pbf
{
    class ParticleSystem : public atlas::utils::Geometry
    {
    public:
        ParticleSystem();

        void setPosition(atlas::math::Point const& pos);

        void updateGeometry(atlas::core::Time<> const& t) override;
        void renderGeometry(atlas::math::Matrix4 const& projection,
            atlas::math::Matrix4 const& view) override;

        
        void resetGeometry() override;
        
        
        // Start
        void init();
        //update grid
        void updateHashPositions();
        bool isValidCell(glm::vec3);

        //step 1
        void addForces(atlas::core::Time<> const& t);
        //step 2 find neighbor
        //Function to return a list of all the neighbors within the specified distance
        // Stored as a pair of index and vector to the neighboring particle
        void findNeighbors(int index);
        
        
        //inside solver=====
        void solveConstraints(atlas::core::Time<> const& t);
        
        void calLambda(int index);
        float poly6Kernel(const glm::vec3 &r, const float h);
        glm::vec3 spikyKernelGradient(const glm::vec3 &r, const float h);
        float gradientConstraints(int index);
        float densityConstraints(int index);

        float densityEstimation(int index);
        
        glm::vec3 calDeltaPi(int index);
        float calScorr(const glm::vec3 &r, const float h);
        void particleCollision(int index,atlas::core::Time<> const& t);
        
        void updatePredPosition(int index);

        void updateVelocity(int index, atlas::core::Time<> const& t);
        
        void updatePosition(int index);

        void applyViscosity(int index);
        float viscosityKernel(const glm::vec3 &r, const float h);
        
        void applyVorticity(int index);
        //end solver=====

        
        
        float bound = 5.f;

        
        //util
        Particle getParticle(int index);
        void addParticles();

    private:
        Grid uniGrid;
        atlas::gl::Buffer mVertexBuffer;
        atlas::gl::Buffer mIndexBuffer;
        atlas::gl::VertexArrayObject mVao;
        GLsizei mIndexCount;
        
        
        
        
        //
        std::vector <Particle> particles;   //List of all the particles in the system
        std::vector <glm::vec3> particlesPos; //List of all the particles' positions
        
        const float poly6KernelConst = 315.0 / (64 * PI);
        const float spikyKernelConst = 45.0 / (PI);
        const float viscosityKernelConst = 15.0 / (2 * PI);
        
        
//        density = 20000 h = 1.
        const float density0 = 10000.0f; //1000kg/m3
        const float smoothingRadius = .5f;
        const int solverIterations = 4;
        
        //Tensile Instability
        const float deltaQ = 0.2 * smoothingRadius;
        const float k = 0.1;
        const int n = 4;

        // for viscosity
        const float c = 0.01;
        const float relaxation = 2.0;

        
        
//        float timeStep = 0.016f;
        glm::vec3 forces = glm::vec3(0.f,-10.f,0.f);

        
//        float cellSize;
        std::map<int, std::vector<int>> hashGrid;
        glm::vec3 gridDim;
        
        

    };
}
