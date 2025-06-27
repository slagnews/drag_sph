#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>

#include "../include/sph.h"
#include "../include/shaders.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void processInput(GLFWwindow *window)
{
    if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

std::vector<float> GetParticleVertices(const std::vector<double>& pos_x, const std::vector<double>& pos_y, const std::vector<double>& vel_x, double Lx, double Ly, double v0) {
    std::vector<float> particleVertices;
    for (size_t i=0; i<pos_x.size(); ++i) {
        float x = static_cast<float>(pos_x[i] / Lx * 2.0 - 1.0);
        float y = static_cast<float>(pos_y[i] / Ly * 2.0 - 1.0);
        float color = static_cast<float>(vel_x[i]) / (2*v0);
        color = std::clamp(color, 0.0f, 1.0f);
        particleVertices.push_back(x);
        particleVertices.push_back(y);
        particleVertices.push_back(color);
    }
    return particleVertices;
}

int main()
{
    double res=0.1;
    int no_steps=2000;
    double Lx=16.0;
    double Ly=16.0;
    double rho0=1.0;
    double kappa=1.2;
    double v0=10e-6;
    double dt=500;
    double c_s=50e-6;
    double beta=0.07;
    double central_radius=1;
    double inflow_factor=0.1;
    double outflow_factor=0.1;
    double kinematic_viscosity=1e-6;
    double max_beta=1.5;

    SimulationParams params(res, no_steps, Lx, Ly, rho0, kappa, v0, dt, c_s, 
        beta, central_radius, inflow_factor, outflow_factor, kinematic_viscosity, max_beta);
    ParticleList pl(params);

    pl.init_grid();
    pl.init_mass();
    pl.init_vel();
    pl.init_type();

    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(800, 600, "LearnOpenGL", NULL, NULL);

    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glViewport(0, 0, 800, 600);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // Get shaders
    Shader ourShader("shaders/basic.vert", "shaders/basic.frag");

    // Retrieve particle positions and send to GPU
    std::vector<float> particleVertices = GetParticleVertices(pl.pos_x, pl.pos_y, pl.vel_x, params.Lx, params.Ly, params.v0);
    unsigned int particleVAO, particleVBO;
    glGenVertexArrays(1, &particleVAO);
    glGenBuffers(1, &particleVBO);
    glBindVertexArray(particleVAO);
    glBindBuffer(GL_ARRAY_BUFFER, particleVBO);
    glBufferData(GL_ARRAY_BUFFER, particleVertices.size()*sizeof(float), particleVertices.data(), GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)(2*sizeof(float)));
    glEnableVertexAttribArray(1);


    while (!glfwWindowShouldClose(window))
    {
        processInput(window);

        glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        ourShader.use();

        // Update particle positions
        pl.integrate();
        particleVertices = GetParticleVertices(pl.pos_x, pl.pos_y, pl.vel_x, params.Lx, params.Ly, params.v0);
        glBindBuffer(GL_ARRAY_BUFFER, particleVBO);
        glBufferData(GL_ARRAY_BUFFER, particleVertices.size() * sizeof(float), particleVertices.data(), GL_DYNAMIC_DRAW);

        // Draw Particles
        float ndc_radius = static_cast<float>(params.central_radius / params.Lx * 2.0);
        ourShader.setFloat("u_radius", ndc_radius);
        glBindVertexArray(particleVAO);
        glPointSize(4.0f);
        glDrawArrays(GL_POINTS, 0, particleVertices.size()/3);
 
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glDeleteVertexArrays(1, &particleVAO);
    glDeleteBuffers(1, &particleVBO);

    glfwTerminate();

    return 0;

}