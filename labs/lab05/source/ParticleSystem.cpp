#include "ParticleSystem.hpp"
#include "Paths.hpp"
#include "LayoutLocations.glsl"

#include <atlas/utils/Mesh.hpp>
#include <atlas/core/STB.hpp>
// Debug
#include <string>
#include <glm/gtx/string_cast.hpp>
#include <iostream>

#include <tbb/parallel_for.h>
#include <tbb/tbb.h>

namespace pbf
{
    ParticleSystem::ParticleSystem() :
        uniGrid(),
        mVertexBuffer(GL_ARRAY_BUFFER),
        mIndexBuffer(GL_ELEMENT_ARRAY_BUFFER)
    {        
        namespace gl = atlas::gl;
        using atlas::math::Point;
        using atlas::math::Vector;
//        printf("%f",particlesPos.size());
        mVao.bindVertexArray();
        mVertexBuffer.bindBuffer();
        mVertexBuffer.bufferData(gl::size<glm::vec3>(particlesPos.size()), particlesPos.data(), GL_DYNAMIC_DRAW);
        mVertexBuffer.vertexAttribPointer(VERTICES_LAYOUT_LOCATION, 3, GL_FLOAT,
                                          GL_FALSE, 0, gl::bufferOffset<float>(0));
        mVao.enableVertexAttribArray(VERTICES_LAYOUT_LOCATION);
        
        mVertexBuffer.unBindBuffer();
        mVao.unBindVertexArray();
        
        std::vector<gl::ShaderUnit> shaders
        {
            {std::string(pbf::ShaderDirectory) + "Spline.vs.glsl", GL_VERTEX_SHADER},
            {std::string(pbf::ShaderDirectory) + "Spline.fs.glsl", GL_FRAGMENT_SHADER}
        };
        
        mShaders.emplace_back(shaders);
        mShaders[0].setShaderIncludeDir(pbf::ShaderDirectory);
        mShaders[0].compileShaders();
        mShaders[0].linkShaders();
        
        auto var = mShaders[0].getUniformVariable("model");
        mUniforms.insert(UniformKey("model", var));
        var = mShaders[0].getUniformVariable("projection");
        mUniforms.insert(UniformKey("projection", var));
        var = mShaders[0].getUniformVariable("view");
        mUniforms.insert(UniformKey("view", var));
        var = mShaders[0].getUniformVariable("colour");
        mUniforms.insert(UniformKey("colour", var));
        
        mShaders[0].disableShaders();
        init();

    }

    
    //start
    void ParticleSystem::init()
    {
        //create particle cube
        glm::vec3 pos(0);
        glm::vec3 v(0);
        // init particles states
        
        float ratio = 5.f/(float)PARTICLE_DIM;
        for( int i = 0 ; i < PARTICLE_DIM; i++)
        {
            for( int j = 0 ; j < PARTICLE_DIM; j++)
            {
                for( int k = 0 ; k < PARTICLE_DIM; k++)
                {
                    float rd1=(rand()%1000)/(10000.0f);// in the range 0 to 999/10000
                    float rd2=(rand()%1000)/(10000.0f);
                    float rd3=(rand()%1000)/(10000.0f);

//                    pos = glm::vec3(i,j,k)* ratio;
                    pos = glm::vec3(i/(rd1+1),j/(rd2+1),k/(rd3+1)) * ratio + glm::vec3(0,.5,0);
                    particles.push_back(Particle(pos, v));
                    particlesPos.push_back(pos);
                    //printf("POSITION %d %d %d \n",i,j,k);
                    //std::cout<<glm::to_string(pos)<<std::endl;
                    
/**      NOTE:   if we use Particle p(pos, v); c++ error
                call to implicity-deleted copy consructor of Particle
                    *p is a pointer and cannot be pushed back.
                   Particle *p = new Particle(pos,v);
                    particles.push_back(&p);**/
                }
            }
        }
        
        //init grid
//        uniGrid.setCellSize(1.5f);
        uniGrid.setCellSize(smoothingRadius);
        uniGrid.setDimension(glm::vec3(40.f));//size of whole cube
        uniGrid.setCenter(glm::vec3(0.f,0.f,0.f));
        uniGrid.setMinExtent(glm::vec3(-20.f));//half dimension
        uniGrid.setMaxExtent(glm::vec3(20.f));
        gridDim = (uniGrid.getMaxExtent() - uniGrid.getMinExtent()) / uniGrid.getCellSize(); // 3*3*3
        for(int i = 0; i < particles.size(); i++)
        {
            if(particles.at(i).getPosition().y>5)
            {
                printf("\n INIT=====\n");
                std::cout<<glm::to_string(particles.at(i).getPosition())<<std::endl;
                printf("%d",i);
                printf("=====\n");
                
                
            }
        }

        
    }
    void ParticleSystem::updateHashPositions()
    {
        hashGrid.clear();
//        tbb::parallel_for(static_cast<std::size_t>(0), particles.size(), [=](const int i)
//        {
//            glm::vec3 hashPosition;
//            glm::vec3 hashPos = glm::floor(
//                                           (particles.at(i).getPredictedPosition() - uniGrid.getMinExtent())/ uniGrid.getCellSize());
//            //update hash positions
//            particles.at(i).setHashPosition(hashPos);
//            std::vector<glm::vec3> neighborCells;
//            
//            neighborCells.push_back(hashPos);
//            for (int xIndex = -1; xIndex <= 1; xIndex++)
//            {
//                for(int yIndex = -1; yIndex <= 1; yIndex++)
//                {
//                    for(int zIndex = -1; zIndex <= 1; zIndex++)
//                    {
//                        neighborCells.push_back(glm::vec3(xIndex,yIndex,zIndex) + hashPos);
//                    }
//                }
//            }
//            
//            for(int j = 0; j < neighborCells.size(); j++)
//            {
//                if(isValidCell(neighborCells.at(j)))
//                {
//                    int cellIndex = neighborCells.at(j).x + gridDim.x * (neighborCells.at(j).y + gridDim.y * neighborCells.at(j).z);
//                    hashGrid[cellIndex].push_back(i);
//                }
//            }
//        });
//        CPU based -by yaoyao
//
        for (int i = 0; i < particles.size(); i++)
        {
            glm::vec3 hashPosition;
            glm::vec3 hashPos = glm::floor((particles.at(i).getPredictedPosition() - uniGrid.getMinExtent())/ uniGrid.getCellSize());
            //update hash positions
            particles.at(i).setHashPosition(hashPos);
            std::vector<glm::vec3> neighborCells;

            neighborCells.push_back(hashPos);
            for (int xIndex = -1; xIndex <= 1; xIndex++)
            {
                for(int yIndex = -1; yIndex <= 1; yIndex++)
                {
                    for(int zIndex = -1; zIndex <= 1; zIndex++)
                    {
                        neighborCells.push_back(glm::vec3(xIndex,yIndex,zIndex) + hashPos);
                    }
                }
            }

            for(int j = 0; j < neighborCells.size(); j++)
            {
                if(isValidCell(neighborCells.at(j)))
                {
                    int cellIndex = neighborCells.at(j).x + gridDim.x * (neighborCells.at(j).y + gridDim.y * neighborCells.at(j).z);
                    hashGrid[cellIndex].push_back(i);
                }
            }
        }
    
    }
    
    bool ParticleSystem::isValidCell(glm::vec3 cellForCheck)
    {
        if(cellForCheck.x >= 0 && cellForCheck.x < gridDim.x)
        {
            if(cellForCheck.y >= 0 && cellForCheck.y < gridDim.y)
            {
                if(cellForCheck.z >= 0 && cellForCheck.z < gridDim.z)
                {
                    return true;
                }
            }
        }
        
        return false;
    }
    
    // Apply external force
    // predicted_position += delta_t*v_i;
    void ParticleSystem::addForces(atlas::core::Time<> const& t)
    {
        tbb::parallel_for(static_cast<std::size_t>(0), particles.size(), [=](const int i)
        {
            Particle & currentParticle = particles.at(i);
            //注意！是forces还是particle.forces？？？
            currentParticle.setVelocity(currentParticle.getVelocity() + t.deltaTime * forces);
            glm::vec3 currPosition = currentParticle.getPosition();
            glm::vec3 predPosition = currPosition + t.deltaTime * currentParticle.getVelocity();
            currentParticle.setPredictedPosition(predPosition);
        });
//        for(int i = 0; i < particles.size(); i++)
////        vector::at is bound-checked and signals if the requested position is out of range by throwing an out_of_range exception."
//        {
//            Particle & currentParticle = particles.at(i);
//            //注意！是forces还是particle.forces？？？
//            currentParticle.setVelocity(currentParticle.getVelocity() + t.deltaTime * forces);
//            glm::vec3 currPosition = currentParticle.getPosition();
//            glm::vec3 predPosition = currPosition + t.deltaTime * currentParticle.getVelocity();
//            currentParticle.setPredictedPosition(predPosition);
//          
//        }
        
    
    }
    
    void ParticleSystem::findNeighbors(int index)
    {
        particles.at(index).clearNeighbors();

        glm::vec3 hashPosition = particles.at(index).getHashPosition();
        int cellIndex = hashPosition.x +
                        gridDim.x * (hashPosition.y + gridDim.y * hashPosition.z);
        //check neighbors in same cell
        int neighborIndex;
        if (hashGrid.find(cellIndex) != hashGrid.end())
        {
            for( int i = 0; i < hashGrid.at(cellIndex).size(); i++)
            {
                neighborIndex = hashGrid.at(cellIndex).at(i);
                if( neighborIndex != index) //if is neighbor not itself
                {
                    float distance = glm::length(particles.at(index).getPredictedPosition()
                                           - particles[neighborIndex].getPredictedPosition());
                    if( distance < smoothingRadius )
                    {
                        particles.at(index).addNeighborIndex(neighborIndex);

                    }
                }
            }
        }
//        else printf("aaaaaaa!!!!! %d %d",cellIndex,index);
    }

    void ParticleSystem::solveConstraints(atlas::core::Time<> const& t)
    {
        tbb::parallel_for(static_cast<std::size_t>(0), particles.size(), [=](const int i)
        {
            calLambda(i);
//            glm::vec3 deltaPi = calDeltaPi(i);
//            particles.at(i).setDeltaPi(deltaPi);
//            particleCollision(i,t);
//            updatePredPosition(i);
//            //clear deltapi
//            particles.at(i).setDeltaPi(glm::vec3(0));
            
        });
        
        tbb::parallel_for(static_cast<std::size_t>(0), particles.size(), [=](const int i)
        {
            glm::vec3 deltaPi = calDeltaPi(i);
            particles.at(i).setDeltaPi(deltaPi);
        });
        
        tbb::parallel_for(static_cast<std::size_t>(0), particles.size(), [=](const int i)
        {
            particleCollision(i,t);
            updatePredPosition(i);
            //clear deltapi
            particles.at(i).setDeltaPi(glm::vec3(0));

        });
        
        tbb::parallel_for(static_cast<std::size_t>(0), particles.size(), [=](const int i)
        {
            updatePredPosition(i);
        });

//        for(int i = 0; i < particles.size(); i++)
//        {
//            calLambda(i);
//            
//            Particle & currentParticle = particles.at(i);
//            glm::vec3 deltaPi = calDeltaPi(i);
//            currentParticle.setDeltaPi(deltaPi);
//            
//            particleCollision(i,t);
//            
//            //update position
//            updatePredPosition(i);
//        
//        }
    
    }

    
    
    void ParticleSystem::calLambda(int index)
    {
        float gradientCstr = gradientConstraints(index);
        // estimate/update density
        float densityI = densityEstimation(index);
        float densityCstr = densityConstraints(index);
        //get lambda i
        float lambdaI = -1.0f * ( densityCstr / ( gradientCstr + relaxation ));

        // update particle's density and lambda
        particles.at(index).setDensity(densityI);
        particles.at(index).setLambda(lambdaI);

    }
    
    float ParticleSystem::poly6Kernel(const glm::vec3 &r, const float h)
    {
        if ( glm::length(r) >= 0.f && glm::length(r) <= h)
        {
            float temp = h * h - glm::length(r) * glm::length(r);
            return (poly6KernelConst / (float)pow(h,9)) * temp * temp * temp;
        }
        else
        {
            return 0.f;
        }
    }
    
    glm::vec3 ParticleSystem::spikyKernelGradient(const glm::vec3 &r, const float h)
    {
        // change from >= 0 to > 0
        if ( glm::length(r) > 0.f && glm::length(r) <= h)
        {
            float temp = (h - glm::length(r)) * (h - glm::length(r));
            return (spikyKernelConst /(float)pow(h,6)) * temp * glm::normalize(r);
        }
        else
        {
            return glm::vec3(0.f,0.f,0.f);
        }
    }

    float ParticleSystem::gradientConstraints(int i)
    {
        std::vector<int> neighbors = particles.at(i).getNeighborIndices();
        int neighborsNum = (int)(neighbors.size());
        int numNeighbors = static_cast<int>(neighbors.size());
        float sumGradient = 0;
        //把所有的position改成了getPredictedPosition
        // get gradient constraints for neighbors
        
        
        for (int j = 0; j < neighbors.size(); j++)
        {
            glm::vec3 r = particles.at(i).getPredictedPosition() - particles.at(neighbors.at(j)).getPredictedPosition();
            glm::vec3 gradientNeighbor = - 1.0f / density0 * spikyKernelGradient(r,smoothingRadius);
            float gradientNeighborLength = glm::length(gradientNeighbor);
            sumGradient += gradientNeighborLength * gradientNeighborLength;
        }
        
        // get gradient constraints for particle itself
        glm::vec3 gradientAtParticle(0,0,0);
        for(int j = 0; j < neighbors.size(); j++)
        {
            glm::vec3 sumGradient(0,0,0);
            glm::vec3 r = particles.at(i).getPredictedPosition() - particles.at(neighbors.at(j)).getPredictedPosition();
            gradientAtParticle += spikyKernelGradient(r,smoothingRadius);
        }
        gradientAtParticle = 1.0f / density0 * gradientAtParticle;
        float gradientAtParticleLength = glm::length(gradientAtParticle);
        sumGradient += gradientAtParticleLength * gradientAtParticleLength;

        return sumGradient;
    }
    
    float ParticleSystem::densityConstraints(int index)
    {
        return  particles.at(index).getDensity() / density0  - 1.0f;
    }

    float ParticleSystem::densityEstimation(int i)
    {
        float density = 0;
        std::vector<int> neighbors = particles.at(i).getNeighborIndices();
        
//        Kernel function implementation
//          Using poly6 kernel function
//        ro += kv * mass;
//        TODO
//          Currently mass is 1 for all particles
//          If mass is changed, change the for loop to multiply by mass
        for( int j = 0; j < neighbors.size(); j++ )
        {
            glm::vec3 r = particles.at(i).getPredictedPosition() - particles.at(neighbors.at(j)).getPredictedPosition();
            density += poly6Kernel(r, smoothingRadius);
        }
        return density;
    }

    
    
    glm::vec3 ParticleSystem::calDeltaPi(int i)
    {
        std::vector<int> neighbors = particles.at(i).getNeighborIndices();
        float lambdaI = particles.at(i).getLambda();
        glm::vec3 deltaPi(0);
        for( int j = 0; j < neighbors.size(); j++ )
        {
            float lambdaJ = particles.at(neighbors.at(j)).getLambda();
            glm::vec3 r = particles.at(i).getPredictedPosition()
                        - particles.at(neighbors.at(j)).getPredictedPosition();

            float sCorr = calScorr(r, smoothingRadius);
            deltaPi += (lambdaI + lambdaJ +sCorr) * spikyKernelGradient(r, smoothingRadius);
        }
        deltaPi /= density0;
        return deltaPi;
    }
    
    float ParticleSystem::calScorr(const glm::vec3 &r, const float h)
    {
        float sCorr = 0.f;
        float Wup = poly6Kernel(r, h);
        float Wdown = poly6Kernel(glm::vec3(deltaQ, 0, 0), h);
//        Wdown = (Wdown == 0.0f) ? 0.00000001 : Wdown;
//        sCorr = -k * pow( Wup / Wdown , n );
        sCorr = ( Wdown < 0.000001f) ? 0.0f : (-k * (float)pow(Wup / Wdown,n));
        return sCorr;
    }

    
    float ParticleSystem::viscosityKernel(const glm::vec3 &r, const float h)
    
    {
        float rLength = glm::length(r);
        if ( rLength >= 0.f && rLength <= h)
        {
            float temp = - (float)pow(rLength / h,3) / 2
                        +  (float)pow(rLength / h,2)
                        + h / (2 * rLength) - 1.0;
            return (viscosityKernelConst / (float)pow(h,3)) * temp;
        }
        else
        {
            return 0.f;
        }
    }

    void ParticleSystem::applyViscosity(int i)
    {
        std::vector<int> neighbors = particles.at(i).getNeighborIndices();

//      sum(j)vij*W
        glm::vec3 sum(0);
//        float c = 1.5; //or 1/density0
        for( int j = 0; j < neighbors.size(); j++ )
        {
            
            sum += (particles.at(neighbors[j]).getVelocity() - particles.at(i).getVelocity()) * viscosityKernel((particles.at(i).getPredictedPosition() - particles.at(neighbors.at(j)).getPredictedPosition()), smoothingRadius);

        }
        glm::vec3 velocity = particles.at(i).getVelocity() + c * sum;
        particles.at(i).setVelocity(velocity);
        
    }

    void ParticleSystem::applyVorticity(int i)
    {
        std::vector<int> neighbors = particles.at(i).getNeighborIndices();
        glm::vec3 omega(0);

        for( int j = 0; j < neighbors.size(); j++ )
        {
            glm::vec3 r = particles.at(i).getPredictedPosition() - particles.at(neighbors.at(j)).getPredictedPosition();
            glm::vec3 vij = particles.at(neighbors.at(j)).getVelocity() - particles.at(i).getVelocity();
            omega += glm::cross(vij, spikyKernelGradient(r, smoothingRadius));
//            calculate a corrective force
            particles.at(i).setVorticity(omega);
        }
        glm::vec3 eta(0);
        for( int j = 0; j < neighbors.size(); j++ )
        {
//            glm::vec3 omegaR = particles.at(i).getVorticity() - particles.at(neighbors.at(j)).getVorticity();
            glm::vec3 r = particles.at(i).getPredictedPosition() - particles.at(neighbors.at(j)).getPredictedPosition();
            eta += spikyKernelGradient(r, smoothingRadius)
                            * glm::length(particles.at(i).getVorticity());
        }
        glm::vec3 N = glm::normalize(eta);
        glm::vec3 force = 0.1f * glm::cross(N, particles.at(i).getVorticity());
        

        
        

    }

    
    // box collision -tbc
    void ParticleSystem::particleCollision(int index , atlas::core::Time<> const& t)
    {
        Particle & ptc = particles.at(index);
        glm::vec3 positionNew = ptc.getPredictedPosition() + ptc.getDeltaPi();
        glm::vec3 positionOld = ptc.getPosition();

        float cr = 0.9;
//        float cf = 0.8;
        float cf = 0.1;

        
        float radius = 0.05f;
        float dampingFactor = 0.02f;
        float eps = 0.0001;
//        x- = x+
        std::vector<glm::vec3> normList = {glm::vec3(0, 1.f, 0),
                                        glm::vec3(-1.f, 0, 0), glm::vec3(1.f, 0 ,0),
                                        glm::vec3(0, 0, -1.f), glm::vec3(0, 0, 1.f)};
        
        // y = 0
        if(positionNew.y + radius < eps)
        {
            cr = 0.5;
            cf = 0.5;
            glm::vec3 norm = normList[0];
            glm::vec3 v = ptc.getVelocity();
            glm::vec3 vN = glm::dot(v, norm) * norm;
            glm::vec3 vT = v - vN;
            glm::vec3 vBounce = -cr * vN;
            glm::vec3 vFriction = (1 - cf) * vT;
            glm::vec3 vNew = vBounce + vFriction;
            
            ptc.setVelocity(vNew);
//            ptc.setVelocity(glm::clamp(vNew,glm::vec3(-0.01,-0.01,-0.01),glm::vec3()));

//            ptc.setPredictedPosition();
        }
        
//        x axis
        if(positionNew.x + radius > bound - eps)
        {
            glm::vec3 norm = normList[1];
            glm::vec3 v = ptc.getVelocity();
            glm::vec3 vN = glm::dot(v, norm) * norm;
            glm::vec3 vT = v - vN;
            glm::vec3 vBounce = -cr * vN;
            glm::vec3 vFriction = (1 - cf) * vT;
            glm::vec3 vNew = vBounce + vFriction;
            ptc.setVelocity(vNew);
        }
//        if(positionNew.x - radius < -bound + eps)
        if(positionNew.x + radius < -0.5-eps)

        {
            glm::vec3 norm = normList[2];
            glm::vec3 v = ptc.getVelocity();
            glm::vec3 vN = glm::dot(v, norm) * norm;
            glm::vec3 vT = v - vN;
            glm::vec3 vBounce = -cr * vN;
            glm::vec3 vFriction = (1 - cf) * vT;
            glm::vec3 vNew = vBounce + vFriction;
            ptc.setVelocity(vNew);
        }
        if(positionNew.z + radius > bound - eps)
        {
            glm::vec3 norm = normList[3];
            glm::vec3 v = ptc.getVelocity();
            glm::vec3 vN = glm::dot(v, norm) * norm;
            glm::vec3 vT = v - vN;
            glm::vec3 vBounce = -cr * vN;
            glm::vec3 vFriction = (1 - cf) * vT;
            glm::vec3 vNew = vBounce + vFriction;
            ptc.setVelocity(vNew);
        }
//        if(positionNew.z - radius < -bound + eps)
        if(positionNew.z + radius < -0.5-eps)

        {
            glm::vec3 norm = normList[4];
            glm::vec3 v = ptc.getVelocity();
            glm::vec3 vN = glm::dot(v, norm) * norm;
            glm::vec3 vT = v - vN;
            glm::vec3 vBounce = -cr * vN;
            glm::vec3 vFriction = (1 - cf) * vT;
            glm::vec3 vNew = vBounce + vFriction;
            ptc.setVelocity(vNew);
        }
        ptc.setPredictedPosition(ptc.getPosition() + t.deltaTime * ptc.getVelocity());
//
//
////        condition for rest
        if (glm::length(ptc.getVelocity()) < EPSILON)
        {
            printf("\n hahaha %d\n",index);
//
        }
        if (ptc.getPredictedPosition().y > 3){
                        printf("\n cao %d\n",index);

                std::cout<<glm::to_string(ptc.getVelocity())<<std::endl;
                std::cout<<glm::to_string(ptc.getPredictedPosition())<<std::endl;
        }

        glm::vec3 newPos = glm::clamp(particles.at(index).getPredictedPosition(),
                                      glm::vec3(-0.2f, -0.0001f, -0.2f),
                                      glm::vec3(bound+0.2, 80, bound+0.2));
        ptc.setPredictedPosition(newPos);

//        py = jt;
//        vec3 normal = vec3(0,1,0);
//        vec3 reflectedDir = particle.velocity - vec3(2.0f*(normal*(Dot(particle.velocity,normal))));
//        particle.velocity[1] = reflectedDir[1]*collision_restitution;
//        return true;

//        for (int c = 0; c < 3; ++c) {
//            if (positionNew[c] <= -bound || (pos[c] >= world_sz_dim[c] - kFloatEpsilon)) {
//                vel[c] = 0.0f;
//                pos[c] =
//                std::max(0.0f, std::min(world_sz_dim[c] - kFloatEpsilon, pos[c]));
//            }
//        }
//
//        ptc.set_position(pos);
//        ptc.set_velocity(vel);
//        
//        float bound = 5.f;
//        Particle & currParticle = particles.at(index);
//
//        glm::vec3 particlePosition = currParticle.getPredictedPosition() + currParticle.getDeltaPi();
////        float radius = 0.05f;
////        float dampingFactor = 0.02f;
////        float eps = 0.000001;
//        if(particlePosition.x - radius < -bound + eps|| particlePosition.x + radius > bound - eps)
//        {
//            
//            currParticle.setVelocity(currParticle.getVelocity() * glm::vec3(-dampingFactor,1,1));
//            currParticle.setPredictedPosition(currParticle.getPosition() + t.deltaTime * currParticle.getVelocity());
//            printf("\n 1 %f----%f %d \n ",uniGrid.getMinExtent().x,uniGrid.getMaxExtent().x,index);
//        }
//        
//        if(particlePosition.y - radius < 0.f + eps)
//        {
//            currParticle.setVelocity(currParticle.getVelocity() * glm::vec3(1,-dampingFactor,1));
//            currParticle.setPredictedPosition(currParticle.getPosition() + t.deltaTime * currParticle.getVelocity());
//            printf("\n 2 %f----%d \n ",uniGrid.getMinExtent().y,index);
//
//        }
//        
//        if(particlePosition.y + radius > uniGrid.getMaxExtent().y - eps)
//        {
//            currParticle.setVelocity(currParticle.getVelocity() * glm::vec3(1,-dampingFactor,1));
//            glm::vec3 pos = currParticle.getPredictedPosition();
//            currParticle.setPredictedPosition(glm::vec3(pos.x, uniGrid.getMaxExtent().y - radius - eps, pos.z));
//            printf("\n 3 %f----%d \n ",uniGrid.getMaxExtent().y,index);
//
//
//        }
//        
//        if(particlePosition.z - radius < -bound + eps || particlePosition.z + radius > bound - eps)
//        {
//            currParticle.setVelocity(currParticle.getVelocity() * glm::vec3(1,1,-dampingFactor));
//            currParticle.setPredictedPosition(currParticle.getPosition() + t.deltaTime * currParticle.getVelocity());
//            printf("\n4 %d\n",index);
//
//        }

        
        
        
        
    }
    
    
    
    
    void ParticleSystem::updatePredPosition(int index)
    {
        glm::vec3 predictedPosition = particles.at(index).getPredictedPosition()
                                        + particles.at(index).getDeltaPi();
        particles.at(index).setPredictedPosition(predictedPosition);
        //if need to clear deltapi?
    
    }
    void ParticleSystem::updateVelocity(int index, atlas::core::Time<> const& t)
    {
        glm::vec3 v = (particles.at(index).getPredictedPosition() - particles.at(index).getPosition()) / t.deltaTime;
        particles.at(index).setVelocity(v);
    
    }

    void ParticleSystem::updatePosition(int index)
    {
        glm::vec3 newPos = glm::clamp(particles.at(index).getPredictedPosition(),
                                      glm::vec3(-0.2f, -0.0001f, -0.2f),
                                      glm::vec3(bound+0.2, 80, bound+0.2));
        
//        glm::vec3 newPos = particles.at(index).getPredictedPosition();
        particles.at(index).setPosition(newPos);
//        particlesPos.push_back(newPos);
//        if (newPos.z < 0.f)
//            printf("woc %d \n",index);
        
    }

    
    void ParticleSystem::setPosition(atlas::math::Point const& pos)
    {
        using atlas::math::Matrix4;
        mModel = glm::translate(Matrix4(1.0f), pos);
        mModel = glm::scale(mModel, glm::vec3(.2f,.2f,.2f));

    }

    void ParticleSystem::updateGeometry(atlas::core::Time<> const& t)
    
    {
        addForces(t);
        updateHashPositions();

//        tbb::parallel_for(0,particles.size(), [=](const int i)
//        {
//            findNeighbors(i);
//        });

        tbb::parallel_for(static_cast<std::size_t>(0), particles.size(), [=](const int i)
        {
            findNeighbors(i);
        });
        
//        for(int i = 0; i < particles.size(); i++)
//        {
//            findNeighbors(i);
//        }

        for(int iter = 0; iter < ITERATIONS; iter++){
            solveConstraints(t);
        }
        
        //update position and velocity

        tbb::parallel_for(static_cast<std::size_t>(0), particles.size(), [=](const int i)
        {
            updateVelocity(i,t);
            applyViscosity(i);
            updatePosition(i);
        });
//        tbb::parallel_for(static_cast<std::size_t>(0), particles.size(), [=](const int i)
//        {
//            applyViscosity(i);
//        });
//        tbb::parallel_for(static_cast<std::size_t>(0), particles.size(), [=](const int i)
//        {
//                        updatePosition(i);
//        });
//        for(int i = 0; i < particles.size(); i++)
//        {
//            updateVelocity(i,t);
//            applyViscosity(i);
//            updatePosition(i);
//            
//        }
        
//        particlesPos.clear();
        printf("==\n");
        std::cout<<glm::to_string(particles.at(0).getPosition())<<std::endl;
        
        particlesPos.clear();
        for(int i = 0; i < particles.size(); i++)
        {
            particlesPos.push_back(particles.at(i).getPosition());
        }
    }
    
    void ParticleSystem::renderGeometry(atlas::math::Matrix4 const& projection,
        atlas::math::Matrix4 const& view)
    {
        namespace math = atlas::math;
        namespace gl = atlas::gl;

        mShaders[0].hotReloadShaders();
        if (!mShaders[0].shaderProgramValid())
        {
            return;
        }

        mShaders[0].enableShaders();

//        mTexture.bindTexture();
        mVao.bindVertexArray();
        mIndexBuffer.bindBuffer();
//        add
        mVertexBuffer.bindBuffer();
        mVertexBuffer.bufferData(gl::size<glm::vec3>(particlesPos.size()), particlesPos.data(), GL_DYNAMIC_DRAW);
        mVertexBuffer.vertexAttribPointer(VERTICES_LAYOUT_LOCATION, 3, GL_FLOAT,
                                          GL_FALSE, 0, gl::bufferOffset<float>(0));
        auto model = mModel;

        glUniformMatrix4fv(mUniforms["model"], 1, GL_FALSE, &mModel[0][0]);
        glUniformMatrix4fv(mUniforms["projection"], 1, GL_FALSE,
            &projection[0][0]);
        glUniformMatrix4fv(mUniforms["view"], 1, GL_FALSE, &view[0][0]);
        
//        const math::Vector sphereColor{ 1.f, .9f, .1f };
//        glUniform3fv(mUniforms["materialColour"], 1, &sphereColor[0]);
//        glDrawElements(GL_TRIANGLES, mIndexCount, GL_UNSIGNED_INT, 0);
        //draw
        glUniform3f(mUniforms["colour"], 200.f, 200.f, 200.f);
        glPointSize(5.f);
        glDrawArrays(GL_POINTS, 0, GLsizei(particlesPos.size()));

        
        mIndexBuffer.unBindBuffer();
        mVao.unBindVertexArray();
//        mTexture.unBindTexture();
        mShaders[0].disableShaders();
    }

    void ParticleSystem::resetGeometry()
    {
    }
    
    void ParticleSystem::addParticles()
    {
        glm::vec3 pos(0);
        glm::vec3 v(0);
        for( int i = 0 ; i < 3; i++)
        {
            for( int j = 0 ; j < 3; j++)
            {
                for( int k = 0 ; k < 3; k++)
                {
                    float rd1=(rand()%1000)/(10000.0f);// in the range 0 to 999/10000
                    float rd2=(rand()%1000)/(10000.0f);
                    float rd3=(rand()%1000)/(10000.0f);
                    
                    pos = glm::vec3(i/(rd1+1),j/(rd2+1),k/(rd3+1)) + glm::vec3(0,.5,0);
                    particles.push_back(Particle(pos, v));
                    particlesPos.push_back(pos);
                  
                }
            }
        }
    
    }

}
