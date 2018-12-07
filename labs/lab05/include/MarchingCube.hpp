#pragma once

#include <atlas/utils/Geometry.hpp>
#include <atlas/gl/Buffer.hpp>
#include <atlas/gl/VertexArrayObject.hpp>

namespace pbf
{
    class MarchingCube : public atlas::utils::Geometry
    {
    public:
        MarchingCube();
////        Grid(glm::vec3 center, glm::vec3 dim,);
//        
////        getter functions
//        glm::vec3 getDimensions();
//        glm::vec3 getCenter();
//        float getCellSize();
//        glm::vec3 getMinExtent();
//        glm::vec3 getMaxExtent();
//        
////        setter functions
//        void setDimension(glm::vec3);
//        void setCenter(glm::vec3);
//        void setCellSize(float);
//        void setMinExtent(glm::vec3);
//        void setMaxExtent(glm::vec3);
//
////        void setPosition(atlas::math::Point const& pos);
////
////        void updateGeometry(atlas::core::Time<> const& t) override;
////        void renderGeometry(atlas::math::Matrix4 const& projection,
////            atlas::math::Matrix4 const& view) override;
////        void resetGeometry() override;
//
//    private:
//        glm::vec3 dimensions;
//        glm::vec3 center;
//        float cellSize;//grid cell size
//        glm::vec3 minExtent;
//        glm::vec3 maxExtent;
//        atlas::gl::Buffer mVertexBuffer;
//        atlas::gl::Buffer mIndexBuffer;
//        atlas::gl::VertexArrayObject mVao;
//        GLsizei mIndexCount;
    };
}
