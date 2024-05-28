#ifndef CUDA_INITS_CUH
#define CUDA_INITS_CUH

#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

#include <application_settings.hpp>

const unsigned int CUDA_TEXTURE_COUNT = 1;
const unsigned int CUDA_RESOURCE_COUNT = 1;
// cudaGraphicsResource_t cuda_resource_1, cuda_resource_2;
cudaGraphicsResource_t *cuda_resources = new cudaGraphicsResource_t[CUDA_RESOURCE_COUNT]; // cuda resource to write to
cudaResourceDesc *cuda_resource_descs = new cudaResourceDesc[CUDA_RESOURCE_COUNT];
GLuint *cuda_textures = new GLuint[CUDA_TEXTURE_COUNT]; // textures used as buffer between CUDA and OpenGL

void createResourceDesciption()
{
    // initialize resource desciption for surface
    for (unsigned int i = 0; i < CUDA_RESOURCE_COUNT; i++)
    {
        memset(&cuda_resource_descs[i], 0, sizeof(cuda_resource_descs[i]));
        cuda_resource_descs[i].resType = cudaResourceTypeArray;
    }

    // memset(&cuda_resource_desc_1, 0, sizeof(cuda_resource_desc_1));
    // cuda_resource_desc_1.resType = cudaResourceTypeArray;
}

void createTextureBuffers(const unsigned int width, const unsigned int height)
{
    glEnable(GL_TEXTURE_2D);

    // ////////////////////////////////////////////////////////////////////////////////
    // // initialize texture 1
    // glGenTextures(1, &textures[0]);
    // glBindTexture(GL_TEXTURE_2D, texture[0]);

    // // initialize texture with glTexImage2D
    // glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);

    // // texture wrapping
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    // // texture filtering
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

    // // TODO implement second texture as buffer for read (OpenGL) and write (CUDA)

    // glBindTexture(GL_TEXTURE_2D, 0);

    // // register OpenGL textures with CUDA
    // checkCudaErrors(cudaGraphicsGLRegisterImage(&cuda_resource_1, texture_1, GL_TEXTURE_2D, cudaGraphicsRegisterFlagsWriteDiscard));

    // createResourceDesciption();

    for (unsigned int i = 0; i < CUDA_TEXTURE_COUNT; i++)
    {
        glGenTextures(1, &cuda_textures[i]);
        glBindTexture(GL_TEXTURE_2D, cuda_textures[i]);

        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
        // texture wrapping
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        // texture filtering
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

        glBindTexture(GL_TEXTURE_2D, 0);
    }

    // register OpenGL textures with CUDA
    for (unsigned int i = 0; i < CUDA_RESOURCE_COUNT; i++)
    {
        checkCudaErrors(cudaGraphicsGLRegisterImage(&cuda_resources[i], cuda_textures[i], GL_TEXTURE_2D, cudaGraphicsRegisterFlagsWriteDiscard));
    }

    // checkCudaErrors(cudaGraphicsGLRegisterImage(&cuda_resources[0], cuda_textures[0], GL_TEXTURE_2D, cudaGraphicsRegisterFlagsWriteDiscard));

    createResourceDesciption();
}

void reset_textures()
{
    for (unsigned int i = 0; i < CUDA_RESOURCE_COUNT; i++)
    {
        cudaGraphicsUnregisterResource(cuda_resources[i]);
    }
    for (unsigned int i = 0; i < CUDA_TEXTURE_COUNT; i++)
    {
        glDeleteTextures(1, &cuda_textures[i]);
    }
}
void cuda_cleanup()
{
    for (unsigned int i = 0; i < CUDA_RESOURCE_COUNT; i++)
    {
        checkCudaErrors(cudaGraphicsUnregisterResource(cuda_resources[i]));
    }
    for (unsigned int i = 0; i < CUDA_TEXTURE_COUNT; i++)
    {
        glDeleteTextures(1, &cuda_textures[i]);
    }
    delete cuda_resources;
    delete cuda_resource_descs;
    delete cuda_textures;
}

#endif