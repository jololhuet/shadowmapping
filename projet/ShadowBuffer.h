#pragma once
#include "icg_common.h"

#ifdef WITH_OPENCV
    #include "opencv2/core/core.hpp"       ///< cv::Mat
    #include "opencv2/highgui/highgui.hpp" ///< cv::imShow
    #include "opencv2/contrib/contrib.hpp" ///< cvtConvert
    #include "opencv2/imgproc/types_c.h"   ///< CV_BGRA2RGBA
#endif

class ShadowBuffer{
protected:
    bool _init;
    int _width;
    int _height;
    GLuint _fbo;
    GLuint _depth_tex;
    GLuint _depthVSM_tex;
    GLint _previous_viewport[4];

public:
    ShadowBuffer(int image_width = 1024, int image_height = 1024){
        this->_width = image_width;
        this->_height = image_height;
    }

    ///--- Warning: overrides viewport!!
    void bind() {
        glGetIntegerv(GL_VIEWPORT, _previous_viewport);  // Store the previous viewport
        glBindFramebuffer(GL_FRAMEBUFFER, _fbo);
        glViewport(0, 0, _width, _height);
        GLenum DrawBuffers[1] = {GL_COLOR_ATTACHMENT0};
        glDrawBuffers(1, DrawBuffers); // "1" is the size of DrawBuffers
    }

    void unbind() {
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glViewport(_previous_viewport[0], _previous_viewport[1], _previous_viewport[2], _previous_viewport[3]);  // Reset the viewport
    }

    void setSize (int image_width, int image_height) {
        this->_width = image_width;
        this->_height = image_height;
    }

    vec2 init(GLint i) {
        ///--- Create (depth, depth²) texture
        {
            glActiveTexture(GL_TEXTURE0+20+i);
            glGenTextures(1, &_depthVSM_tex);
            glBindTexture(GL_TEXTURE_2D, _depthVSM_tex);

            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, _width, _height, 0, GL_RGBA, GL_FLOAT, NULL);

            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

            //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glGenerateMipmap(GL_TEXTURE_2D);
        }

        ///--- Create Depth texture
        {
            glActiveTexture(GL_TEXTURE0+10+i);
            glGenTextures(1, &_depth_tex);
            glBindTexture(GL_TEXTURE_2D, _depth_tex);

            glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, _width, _height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);

            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

            //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        }

        ///--- Tie it all together
        {
            glGenFramebuffers(1, &_fbo);
            glBindFramebuffer(GL_FRAMEBUFFER, _fbo);
            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, _depthVSM_tex, 0);
            glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, _depth_tex, 0);

            if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
                std::cerr << "!!!ERROR: Framebuffer not OK :(" << std::endl;

            glBindFramebuffer(GL_FRAMEBUFFER, 0);
        }

        return vec2(_depth_tex, _depthVSM_tex);
    }

    void cleanup() {
        glBindFramebuffer(GL_FRAMEBUFFER, 0 /*UNBIND*/);
        glDeleteFramebuffers(1, &_fbo);
        glDeleteTextures(1, &_depthVSM_tex);
        glDeleteTextures(1, &_depth_tex);
    }
};
