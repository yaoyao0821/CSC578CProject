#include "UniformGrid.hpp"
#include "Paths.hpp"
#include "LayoutLocations.glsl"

#include <atlas/utils/Mesh.hpp>
#include <atlas/core/STB.hpp>
//#include <glm/gtx/string_cast.hpp>
//#include <iostream>
//#include "glm/ext.hpp"
namespace pbf
{
    Grid::Grid()
//        mVertexBuffer(GL_ARRAY_BUFFER),
//        mIndexBuffer(GL_ELEMENT_ARRAY_BUFFER)
    {
        
    }
    

    


    //getter functions
    float Grid::getCellSize()
    {
        return this->cellSize;
    }
    
    glm::vec3 Grid::getDimensions()
    {
        return this->dimensions;
    }
    
    glm::vec3 Grid::getCenter()
    {
        return this->center;
    }
    
    glm::vec3 Grid::getMinExtent()
    {
        return this->minExtent;
    }
    glm::vec3 Grid::getMaxExtent()
    {
        return this->maxExtent;
    }
    
    //setter functions
    void Grid::setCellSize(float size)
    {
        this->cellSize = size;
    }
    void Grid::setDimension(glm::vec3 dim)
    {
        this->dimensions = dim;
    }
    void Grid::setCenter(glm::vec3 c)
    {
        this->center = c;
    }
    void Grid::setMinExtent(glm::vec3 min)
    {
        this->minExtent = min;
    }
    void Grid::setMaxExtent(glm::vec3 max)
    {
        this->maxExtent = max;
    }
    
//        using atlas::utils::Mesh;
//        namespace gl = atlas::gl;
//        namespace math = atlas::math;
//
//        Mesh sphere;
//        std::string path{ DataDirectory };
//        path = path + "cube.obj";
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
//        mModel = glm::scale(math::Matrix4(1.0f), math::Vector(1.0f));
//        mModel = glm::scale(math::Matrix4(1.0f), glm::vec3(0.5f,0.5f,0.5f));
//        std::cout<<glm::to_string(mVelocity)<<std::endl;
        
//    }

//    void Grid::setPosition(atlas::math::Point const& pos)
//    {
//        using atlas::math::Matrix4;
////        using atlas::math::cartesianToPolar;
//        mModel = glm::translate(Matrix4(1.0f), pos);
////        mModel = glm::scale(mModel, glm::vec3(.4f,.4f,.4f));
//    }

//    void Grid::updateGeometry(atlas::core::Time<> const& t)
//    {
//        using atlas::math::Matrix4;
//        using atlas::math::Vector;
//        using atlas::math::Point;
//        if (glm::abs(mPosition.x) < 6.0f)
//        {
//            return;
//        }
//
//        if (glm::abs(mPosition.x) > 50.0f)
//        {
//            mVelocity.x = -1.0f * mVelocity.x;
//        }

        
//    }
    
//    void Grid::renderGeometry(atlas::math::Matrix4 const& projection,
//        atlas::math::Matrix4 const& view)
//    {
//        namespace math = atlas::math;
//
//        mShaders[0].hotReloadShaders();
//        if (!mShaders[0].shaderProgramValid())
//        {
//            return;
//        }
//
//        mShaders[0].enableShaders();
//
//        mVao.bindVertexArray();
//        mIndexBuffer.bindBuffer();
//
//        auto model = mModel;
//
//        glUniformMatrix4fv(mUniforms["model"], 1, GL_FALSE, &mModel[0][0]);
//        glUniformMatrix4fv(mUniforms["projection"], 1, GL_FALSE,
//            &projection[0][0]);
//        glUniformMatrix4fv(mUniforms["view"], 1, GL_FALSE, &view[0][0]);
//        const math::Vector sphereColor{ 1.f, 0.f, 0.f };
//        glUniform3fv(mUniforms["materialColour"], 1, &sphereColor[0]);
//        
//        glDrawElements(GL_TRIANGLES, mIndexCount, GL_UNSIGNED_INT, 0);
//
//        mIndexBuffer.unBindBuffer();
//        mVao.unBindVertexArray();
//        mShaders[0].disableShaders();
//    }
//
//    void Grid::resetGeometry()
//    {
//    }
}
