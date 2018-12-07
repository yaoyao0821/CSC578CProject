#pragma once

#include "UniformGrid.hpp"
#include "ParticleSystem.hpp"

#include <atlas/tools/ModellingScene.hpp>
#include <atlas/utils/FPSCounter.hpp>
#include <atlas/tools/MayaCamera.hpp>


namespace pbf
{
    class FluidsScene : public atlas::tools::ModellingScene
    {
    public:
        FluidsScene();

        void updateScene(double time) override;
        void renderScene() override;
        //camera
        void mousePressEvent(int button, int action, int modifiers,
                             double xPos, double yPos) override;
        void mouseMoveEvent(double xPos, double yPos) override;
    private:
        atlas::utils::FPSCounter mAnimCounter;
        ParticleSystem ps;
        
        bool mPlay;
        
        bool frontView;
        
    };
}

