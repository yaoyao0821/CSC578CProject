#pragma once

#include <atlas/utils/Geometry.hpp>
#include <atlas/gl/Buffer.hpp>
#include <atlas/gl/VertexArrayObject.hpp>

namespace pbf
{
    class Particle : public atlas::utils::Geometry
    {
    public:
        Particle();
        Particle(glm::vec3);    //Construct a particle at a specific position
        Particle(glm::vec3, glm::vec3); //Construct the particle with pos and vel


//        void setPosition(atlas::math::Point const& pos);

        void updateGeometry(atlas::core::Time<> const& t) override;
        void renderGeometry(atlas::math::Matrix4 const& projection,
            atlas::math::Matrix4 const& view) override;

        void resetGeometry() override;
        
        
        //getter funcitons
        float getMass();
        float getRadius();
        glm::vec3 getPosition();
        glm::vec3 getVelocity();
        glm::vec3 getForces();
        glm::vec3 getPredictedPosition();
        float getLambda();
        glm::vec3 getDeltaPi();
        
        glm::ivec3 getHashPosition();
        float getDensity();
        std::vector<int> getNeighborIndices();
        glm::vec3 getVorticity();

        //setter functions
        void setPosition(glm::vec3);
        void setVelocity(glm::vec3);
        void setForces(glm::vec3);
        void setPredictedPosition(glm::vec3);
        void setLambda(float);
        void setDeltaPi(glm::vec3);
        void setHashPosition(glm::ivec3);
        void setDensity(float d);
        void setVorticity(glm::vec3);

        //neighbors functions
        void addNeighborIndex(int index);
        void clearNeighbors();
        

    private:

        atlas::gl::Buffer mVertexBuffer;
        atlas::gl::Buffer mIndexBuffer;
        atlas::gl::VertexArrayObject mVao;
        GLsizei mIndexCount;
        
        //
        float mass = 1.0f;    //Mass of all particles will be same
        float radius = 0.5f;  //Radius of all particles will be same
        float density = 0.0f;
        
        //  Current values
        glm::vec3 position;     //Stores the current position of the particle
        glm::vec3 velocity;     //Stores the current velocity of the particle
        glm::vec3 forces;       //Stores the accumulated forces acting on the particle
        
        //old position
        glm::vec3 oldPosition;

        
        //predicted values
        glm::vec3 predictedPosition;    //Predicted particle position
        glm::vec3 deltaPi;
        
        glm::ivec3 hashPosition; //stores the grid index where the particle predicted position is
        std::vector <int> neighborIndices;    //stores the indices of all the neighbors
        float lambda;
        
        //vorticity
        glm::vec3 vorticity;

    };
}
