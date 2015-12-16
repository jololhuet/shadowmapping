#pragma once
#include "icg_common.h"
#include "attrib_locations.h"

class FrustPoint {
protected:
    GLuint _pid;
    GLuint _vao;
    vec3 _position;

public:
    FrustPoint() {
        _position = vec3(0,0,0);
    }

    FrustPoint(vec3 pos) {
        _position = pos;
    }

    void init () {

        /// Shader Initialisation
        {
            _pid = opengp::load_shaders("point_vshader.glsl", "point_fshader.glsl");
            glBindAttribLocation(_pid, ATTRIB_LOC_vpoint, "vpoint");
            glLinkProgram(_pid);
        }

        glGenVertexArrays(1, &_vao);
        glBindVertexArray(_vao);

        glBindVertexArray(0);
        glUseProgram(0);
        check_error_gl();
    }

    void draw (const mat4 view, const mat4 projection) {
        glUseProgram(_pid);
        glBindVertexArray(_vao);

        ///--- Load position
        GLint pos_id = glGetUniformLocation(_pid, "position");
        glUniform3fv(pos_id, 1, _position.data());

        mat4 VP = projection*view;
        GLint VP_id = glGetUniformLocation(_pid, "vp");
        assert(VP_id >= 0);
        glUniformMatrix4fv(VP_id, 1, GL_FALSE, VP.data());

        glEnable(GL_PROGRAM_POINT_SIZE);
        glDrawArrays(GL_POINTS, 0, 1);
        glDisable(GL_PROGRAM_POINT_SIZE);

        glBindVertexArray(0);
        glUseProgram(0);

        check_error_gl();
    }

    void cleanup() {
        glDeleteVertexArrays(1, &_vao);
    }

    void setPosition (vec3 newPos) {
        _position = newPos;
    }

};


