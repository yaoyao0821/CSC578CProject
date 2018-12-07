#include "FluidsScene.hpp"

#include <atlas/utils/GUI.hpp>
#include <atlas/gl/GL.hpp>

namespace pbf
{
    FluidsScene::FluidsScene() :
//        mSpring(atlas::math::Vector(5.f, 0.0, 0),
//            atlas::math::Vector(0, 0.0f, 5.0f),
//            atlas::math::Vector(0,1,0)),

        mPlay(false),
        frontView(false),
        mAnimCounter(45.0f)
    {}

    void FluidsScene::updateScene(double time)
    {
        using atlas::core::Time;
        
        ModellingScene::updateScene(time);
        if (mPlay && mAnimCounter.isFPS(mTime))
        {
            ps.updateGeometry(mTime);
//            lagrange.updateGeometry(mTime,mPlay);
//            mSpring.setPosition(atlas::math::Vector(0, 0, 0));
//            cube2.setPosition(atlas::math::Vector(0, 5, 0));

        }
        
    }
    //camera
    void FluidsScene::mousePressEvent(int button, int action, int modifiers,
                                      double xPos, double yPos)
    {
        using atlas::tools::MayaMovements;
        atlas::utils::Gui::getInstance().mousePressed(button, action, modifiers);
        
        if (action == GLFW_PRESS)
        {
            atlas::math::Point2 point(xPos, yPos);
            
            if (button == GLFW_MOUSE_BUTTON_LEFT&&
                modifiers == GLFW_MOD_ALT)
            {
                mCamera.setMovement(MayaMovements::Tumble);
                mCamera.mouseDown(point);
            }
            else
            {
                mCamera.mouseUp();
            }
            
        }
    }
    void FluidsScene::mouseMoveEvent(double xPos, double yPos)
    {
        atlas::utils::Gui::getInstance().mouseMoved(xPos, yPos);
        mCamera.mouseMove(atlas::math::Point2(xPos, yPos));
        
    }
    

    void FluidsScene::renderScene()
    {
        using atlas::utils::Gui;
        
        Gui::getInstance().newFrame();
        const float grey = 0.0f / 255.0f;
        glClearColor(grey, grey, grey, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        mProjection = glm::perspective(
            glm::radians(mCamera.getCameraFOV()),
            (float)mWidth / mHeight, 1.0f, 100000000.0f);
        if (frontView)
        {
            using atlas::math::Vector;
            auto cameraPosition = Vector(0.0f,10.0f,40.0f);
            auto cameraTarget = Vector(0.,0.,0.);
            auto mUp = Vector(0, 1, 0);
            mView =  glm::lookAt(cameraPosition, cameraTarget, mUp);
        }
        else
        {
            mView = mCamera.getCameraMatrix();
        }

//        mSpring.renderGeometry(mProjection, mView);
        ps.renderGeometry(mProjection, mView);
        mGrid.renderGeometry(mProjection, mView);
//        lagrange.renderGeometry(mProjection, mView);
//        drawSphere((-1,-1,-1));

        // Global HUD
        ImGui::SetNextWindowSize(ImVec2(350, 140), ImGuiSetCond_FirstUseEver);
        ImGui::Begin("Global HUD");
//        if (ImGui::Button("Reset Camera"))
//        {
//            mCamera.resetCamera();
//        }

        if (mPlay)
        {
            if (ImGui::Button("Pause"))
            {
                mPlay = !mPlay;
            }
        }
        else
        {
            if (ImGui::Button("Play"))
            {
                mPlay = !mPlay;
            }
        }

        if (ImGui::Button("Reset"))
        {
            mPlay = false;
            mTime.currentTime = 0.0f;
            mTime.totalTime = 0.0f;
        }
        if (frontView)
        {
            if (ImGui::Button("Reset Camera"))
            {
                frontView = !frontView;
                mCamera.resetCamera();
            }
        
        }
        else
        {
            if (ImGui::Button("Front View"))
            {
                frontView = !frontView;
            }
        }
        if (ImGui::Button("Add"))
        {
            ps.addParticles();
        }
        

      
        ImGui::Text("Application average %.3f ms/frame (%.1FPS)",
            1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
        ImGui::End();

        ImGui::SetNextWindowSize(ImVec2(300, 150), ImGuiSetCond_FirstUseEver);
        ImGui::Begin("bound");
        ImGui::InputFloat("bound", &ps.bound, 0.1f, 5.0f, 1);

        ImGui::End();
//
//        ImGui::SetNextWindowSize(ImVec2(300, 150), ImGuiSetCond_FirstUseEver);
//        ImGui::Begin("Sphere States (s)");
//        ImGui::InputFloat("s", &lagrange.s, 0.1f, 5.0f, 1);
//        ImGui::InputFloat("sdot", &lagrange.sdot, 0.1f, 5.0f, 1);
//        ImGui::InputFloat("sdotdot", &lagrange.sdotdot, 0.1f, 5.0f, 1);
//        ImGui::End();
//        
//        
//        ImGui::SetNextWindowSize(ImVec2(300, 150), ImGuiSetCond_FirstUseEver);
//        ImGui::Begin("Cube States (x)");
//        ImGui::InputFloat("x", &lagrange.x, 0.1f, 5.0f, 1);
//        ImGui::InputFloat("xdot", &lagrange.xdot, 0.1f, 5.0f, 1);
//        ImGui::InputFloat("xdotdot", &lagrange.xdotdot, 0.1f, 5.0f, 1);
//        ImGui::End();
//        
//        ImGui::SetNextWindowSize(ImVec2(300, 150), ImGuiSetCond_FirstUseEver);
//        ImGui::Begin("Other Settings");
//        ImGui::InputFloat("Mass(cube)", &lagrange.M, 1.0f, 5.0f, 1);
//        ImGui::InputFloat("mass(sphere)", &lagrange.m, 1.0f, 5.0f, 1);
//        ImGui::InputFloat("k(spring)", &lagrange.k, 1.0f, 5.0f, 1);
//        ImGui::InputFloat("l(spring)", &lagrange.l, 1.0f, 10.0f, 1);
//
//        ImGui::End();
        
        ImGui::Render();
    }
    
    
 
    
    
    
    
    
}
