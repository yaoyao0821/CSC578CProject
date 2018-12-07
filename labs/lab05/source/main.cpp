#include "FluidsScene.hpp"

#include <atlas/utils/Application.hpp>
#include <atlas/utils/WindowSettings.hpp>
#include <atlas/gl/ErrorCheck.hpp>
#include <tbb/parallel_for.h>
int main()
{
    using atlas::utils::WindowSettings;
    using atlas::utils::ContextVersion;
    using atlas::utils::Application;
    using atlas::utils::ScenePointer;
    using namespace pbf;

    atlas::gl::setGLErrorSeverity(
        ATLAS_GL_ERROR_SEVERITY_HIGH | ATLAS_GL_ERROR_SEVERITY_MEDIUM);

    WindowSettings settings;
    settings.contextVersion = { 3, 3 };
    settings.isForwardCompat = true;
    settings.isMaximized = true;
//    printf("Ca!");
//    std::vector <int> p;
//    p.push_back(1);
//    p.push_back(2);
//    p.push_back(3);
//
//    tbb::parallel_for(static_cast<std::size_t>(0), p.size(), [=](const int r)
//                      {
//                          printf("%d",p.at(r));
//                      });
//    tbb::parallel_for(0,8,
//                      [=](const int r)
//                      {
//                          printf("%d",r);
//                      });
    Application::getInstance().createWindow(settings);
    Application::getInstance().addScene(ScenePointer(new FluidsScene));
    Application::getInstance().runApplication();

    return 0;

}
