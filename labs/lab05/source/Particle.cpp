#include "Particle.hpp"
#include "Paths.hpp"
#include "LayoutLocations.glsl"

#include <atlas/utils/Mesh.hpp>
#include <atlas/core/STB.hpp>

namespace pbf
{
//    Particle::Particle() :
//        mVertexBuffer(GL_ARRAY_BUFFER),
//        mIndexBuffer(GL_ELEMENT_ARRAY_BUFFER)
//    {
//        using atlas::utils::Mesh;
//        namespace gl = atlas::gl;
//        namespace math = atlas::math;
//
//        Mesh sphere;
//        std::string path{ DataDirectory };
//        path = path + "sphere.obj";
//        Mesh::fromFile(path, sphere);
//
//        mIndexCount = static_cast<GLsizei>(sphere.indices().size());
//
//        std::vector<float> data;
//        for (std::size_t i = 0; i < sphere.vertices().size(); ++i)
//        {
//            data.push_back(sphere.vertices()[i].x);
//            data.push_back(sphere.vertices()[i].y);
//            data.push_back(sphere.vertices()[i].z);
//
//            data.push_back(sphere.normals()[i].x);
//            data.push_back(sphere.normals()[i].y);
//            data.push_back(sphere.normals()[i].z);
//
//            data.push_back(sphere.texCoords()[i].x);
//            data.push_back(sphere.texCoords()[i].y);
//        }
//
//        mVao.bindVertexArray();
//        mVertexBuffer.bindBuffer();
//        mVertexBuffer.bufferData(gl::size<float>(data.size()), data.data(),
//            GL_STATIC_DRAW);
//        mVertexBuffer.vertexAttribPointer(VERTICES_LAYOUT_LOCATION, 3, GL_FLOAT,
//            GL_FALSE, gl::stride<float>(8), gl::bufferOffset<float>(0));
//        mVertexBuffer.vertexAttribPointer(NORMALS_LAYOUT_LOCATION, 3, GL_FLOAT,
//            GL_FALSE, gl::stride<float>(8), gl::bufferOffset<float>(3));
//        mVertexBuffer.vertexAttribPointer(TEXTURES_LAYOUT_LOCATION, 2, GL_FLOAT,
//            GL_FALSE, gl::stride<float>(8), gl::bufferOffset<float>(6));
//
//        mVao.enableVertexAttribArray(VERTICES_LAYOUT_LOCATION);
//        mVao.enableVertexAttribArray(NORMALS_LAYOUT_LOCATION);
//        mVao.enableVertexAttribArray(TEXTURES_LAYOUT_LOCATION);
//
//        mIndexBuffer.bindBuffer();
//        mIndexBuffer.bufferData(gl::size<GLuint>(sphere.indices().size()),
//            sphere.indices().data(), GL_STATIC_DRAW);
//
//        mIndexBuffer.unBindBuffer();
//        mVertexBuffer.unBindBuffer();
//        mVao.unBindVertexArray();
//
//        std::vector<gl::ShaderUnit> shaders
//        {
//            {std::string(ShaderDirectory) + "Ball.vs.glsl", GL_VERTEX_SHADER},
//            {std::string(ShaderDirectory) + "Ball.fs.glsl", GL_FRAGMENT_SHADER}
//        };
//
//        mShaders.emplace_back(shaders);
//        mShaders[0].setShaderIncludeDir(ShaderDirectory);
//        mShaders[0].compileShaders();
//        mShaders[0].linkShaders();
//
//        auto var = mShaders[0].getUniformVariable("model");
//        mUniforms.insert(UniformKey("model", var));
//        var = mShaders[0].getUniformVariable("projection");
//        mUniforms.insert(UniformKey("projection", var));
//        var = mShaders[0].getUniformVariable("view");
//        mUniforms.insert(UniformKey("view", var));
//        var = mShaders[0].getUniformVariable("materialColour");
//        mUniforms.insert(UniformKey("materialColour", var));
//
//        mShaders[0].disableShaders();
//
//
//        
//    }
    
    
    //static variables
    
    //float Particle::mass = 1.0f;
    //float Particle::radius = 1.0f;
    
    //Constructors
    Particle::Particle()
    {
        this->position = glm::vec3(0,0,0);
        this->velocity = glm::vec3(0,0,0);
        this->predictedPosition = glm::vec3(0,0,0);
        this->deltaPi = glm::vec3(0,0,0);
        this->vorticity = glm::vec3(0,0,0);

    }
    
    Particle::Particle(glm::vec3 pos)
    {
        this->position = pos;
        this->velocity = glm::vec3(0,0,0);
        this->predictedPosition = glm::vec3(0,0,0);
        this->deltaPi = glm::vec3(0,0,0);
        this->vorticity = glm::vec3(0,0,0);

    }
    
    Particle::Particle(glm::vec3 pos, glm::vec3 vel)
    {
        this->position = pos;
        this->velocity = vel;
        this->predictedPosition = glm::vec3(0,0,0);
        this->deltaPi = glm::vec3(0,0,0);
        this->vorticity = glm::vec3(0,0,0);

    }
    
    
    
    
    //Getter functions
    float Particle::getMass()
    {
        return this->mass;
    }
    
    float Particle::getRadius()
    {
        return this->radius;
    }
    
    glm::vec3 Particle::getPosition()
    {
        return this->position;
    }
    
    glm::vec3 Particle::getVelocity()
    {
        return this->velocity;
    }
    
    glm::vec3 Particle::getForces()
    {
        return this->forces;
    }
    
    glm::vec3 Particle::getPredictedPosition()
    {
        return this->predictedPosition;
    }
    
    float Particle::getLambda()
    {
        return lambda;
    }
    
    glm::vec3 Particle::getDeltaPi()
    {
        return deltaPi;
    }
    
    float Particle::getDensity(){
        return density;
    }
    
    void Particle::setDensity(float density)
    {
        this->density = density;
    }
    
    std::vector<int> Particle::getNeighborIndices()
    {
        return neighborIndices;
    }
    
    glm::ivec3 Particle::getHashPosition()
    {
        return hashPosition;
    }
    glm::vec3 Particle::getVorticity()
    {
        return vorticity;
    }

    //Setter functions
    void Particle::setPosition(glm::vec3 pos)
    {
        this->position = pos;
    }
    
    void Particle::setVelocity(glm::vec3 vel)
    {
        this->velocity = vel;
    }
    
    void Particle::setForces(glm::vec3 netForce)
    {
        this->forces = netForce;
    }
    
    void Particle::setPredictedPosition(glm::vec3 pos)
    {
        this->predictedPosition = pos;
    }
    
    void Particle::addNeighborIndex(int index)
    {
        neighborIndices.push_back(index);
    }
    
    void Particle::clearNeighbors()
    {
        neighborIndices.clear();
    }
    
    void Particle::setLambda(float l)
    {
        lambda = l;
    }
    
    void Particle::setDeltaPi(glm::vec3 p)
    {
        deltaPi = p;
    }
    
    void Particle::setHashPosition(glm::ivec3 p)
    {
        hashPosition = p;
    }
    
    void Particle::setVorticity(glm::vec3 v)
    {
        this->vorticity = v;
    }

    
    //

//    void Particle::setPosition(atlas::math::Point const& pos)
//    {
//        using atlas::math::Matrix4;
//        mModel = glm::translate(Matrix4(1.0f), pos);
//        mModel = glm::scale(mModel, glm::vec3(.2f,.2f,.2f));
//
//    }

    void Particle::updateGeometry(atlas::core::Time<> const& t)
    
    {}
    
    void Particle::renderGeometry(atlas::math::Matrix4 const& projection,
        atlas::math::Matrix4 const& view)
    {
        namespace math = atlas::math;

        mShaders[0].hotReloadShaders();
        if (!mShaders[0].shaderProgramValid())
        {
            return;
        }

        mShaders[0].enableShaders();

//        mTexture.bindTexture();
        mVao.bindVertexArray();
        mIndexBuffer.bindBuffer();

        auto model = mModel;

        glUniformMatrix4fv(mUniforms["model"], 1, GL_FALSE, &mModel[0][0]);
        glUniformMatrix4fv(mUniforms["projection"], 1, GL_FALSE,
            &projection[0][0]);
        glUniformMatrix4fv(mUniforms["view"], 1, GL_FALSE, &view[0][0]);
        
        const math::Vector sphereColor{ 1.f, .9f, .1f };
        glUniform3fv(mUniforms["materialColour"], 1, &sphereColor[0]);
        
        glDrawElements(GL_TRIANGLES, mIndexCount, GL_UNSIGNED_INT, 0);

        mIndexBuffer.unBindBuffer();
        mVao.unBindVertexArray();
//        mTexture.unBindTexture();
        mShaders[0].disableShaders();
    }

    void Particle::resetGeometry()
    {
    }
    
    
    
    
    
}
