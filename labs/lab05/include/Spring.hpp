#pragma once

#include <atlas/utils/Geometry.hpp>
#include <atlas/gl/Buffer.hpp>
#include <atlas/gl/VertexArrayObject.hpp>
#include <atlas/gl/Texture.hpp>

#include <glm/vec2.hpp>
#include <array>
#include <vector>
#include <cmath>

namespace pbf
{
    class Spring : public atlas::utils::Geometry
    {
        public:
            Spring(atlas::math::Point const& start,atlas::math::Point const& end,atlas::math::Vector colour);
            Spring();


                void updateGeometry(atlas::core::Time<> const& t) override;
                void renderGeometry(atlas::math::Matrix4 const& projection,
                            atlas::math::Matrix4 const& view) override;
                void resetGeometry() override;
                void setPoints(atlas::math::Point const& start,atlas::math::Point const& end);


        private:

            atlas::gl::Buffer mVertexBuffer;
            atlas::gl::Buffer mIndexBuffer;
            atlas::gl::VertexArrayObject mVao;
            std::array<atlas::math::Vector, 2> mPoints;
            atlas::math::Vector colour;


    };
}

