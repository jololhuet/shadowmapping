#pragma once
#include "icg_common.h"
#include "../attrib_locations.h"

class Floor{
protected:
    GLuint _vao; ///< vertex array object
    GLuint _tex; ///< Texture ID
    GLuint _texNormal; ///< normal texture ID

    mat4 model = mat4::Identity();
    mat4 bigModel;

public:
    void init(){
        float s = 4.0;
        mat4 scale;
        scale << s,    0.0f,  0.0f,  0.0f,
                      0.0f,  s,    0.0f,  0.0f,
                      0.0f,  0.0f,  s,    0.0f,
                      0.0f,  0.0f,  0.0f,  1.0f;

        mat4 trans = Eigen::Affine3f(Eigen::Translation3f(vec3(0.0f, 0.0f, 0.0f))).matrix();

        bigModel = trans * scale;
        ///--- Vertex one vertex Array
        glGenVertexArrays(1, &_vao);
        glBindVertexArray(_vao);

        ///--- Vertex coordinates
        {
            const GLfloat vpoint[] = { /*V1*/ -1.0f,  0.0f,  1.0f,
                                       /*V2*/  1.0f,  0.0f,  1.0f,
                                       /*V3*/ -1.0f,  0.0f, -1.0f,
                                       /*V4*/  1.0f,  0.0f, -1.0f };
            ///--- Buffer
            GLuint vbo;
            glGenBuffers(1, &vbo);
            glBindBuffer(GL_ARRAY_BUFFER, vbo);
            glBufferData(GL_ARRAY_BUFFER, sizeof(vpoint), vpoint, GL_STATIC_DRAW);

            ///--- Attribute
            glEnableVertexAttribArray(ATTRIB_LOC_vpoint);
            glVertexAttribPointer(ATTRIB_LOC_vpoint, 3, GL_FLOAT, DONT_NORMALIZE, ZERO_STRIDE, ZERO_BUFFER_OFFSET);
        }

        ///--- Normal coordinates
        {
            const GLfloat vnormal[] = { 0.0f, +1.0f, 0.0f,
                                        0.0f, +1.0f, 0.0f,
                                        0.0f, +1.0f, 0.0f,
                                        0.0f, +1.0f, 0.0f };
            ///--- Buffer
            GLuint vbo;
            glGenBuffers(1, &vbo);
            glBindBuffer(GL_ARRAY_BUFFER, vbo);
            glBufferData(GL_ARRAY_BUFFER, sizeof(vnormal), vnormal, GL_STATIC_DRAW);

            ///--- Attribute
            glEnableVertexAttribArray(ATTRIB_LOC_vnormal);
            glVertexAttribPointer(ATTRIB_LOC_vnormal, 3, GL_FLOAT, DONT_NORMALIZE, ZERO_STRIDE, ZERO_BUFFER_OFFSET);
        }

        ///--- Texture coordinates
        {
            const GLfloat vtexcoord[] = { /*V1*/ 0.0f, 0.0f, 
                                          /*V2*/ 1.0f, 0.0f, 
                                          /*V3*/ 0.0f, 1.0f,
                                          /*V4*/ 1.0f, 1.0f};

            ///--- Buffer
            GLuint vbo;
            glGenBuffers(1, &vbo);
            glBindBuffer(GL_ARRAY_BUFFER, vbo);
            glBufferData(GL_ARRAY_BUFFER, sizeof(vtexcoord), vtexcoord, GL_STATIC_DRAW);

            ///--- Attribute
            glEnableVertexAttribArray(ATTRIB_LOC_vtexcoord);
            glVertexAttribPointer(ATTRIB_LOC_vtexcoord, 2, GL_FLOAT, DONT_NORMALIZE, ZERO_STRIDE, ZERO_BUFFER_OFFSET);
        }

        check_error_gl();

        ///--- Load texture
        glGenTextures(1, &_tex);
        glBindTexture(GL_TEXTURE_2D, _tex);
        check_error_gl();
        if (!glfwLoadTexture2D("planks.tga", 0)) {
          std::cerr << "Couldn't load planks.tga" << std::endl;
        }
        check_error_gl();
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glGenerateMipmap(GL_TEXTURE_2D);

        check_error_gl();

        glGenTextures(1, &_texNormal);
        glBindTexture(GL_TEXTURE_2D, _texNormal);
        if (!glfwLoadTexture2D("planks_normal.tga", 0)) {
          std::cerr << "Couldn't load planks_normal.tga" << std::endl;
        }
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glGenerateMipmap(GL_TEXTURE_2D);
        check_error_gl();
    }

    /**
     * @brief setSize
     * @param size 0 -> small ; 1 (or !0) -> big
     */
    void setSize (int size) {
        if (size == 0) {
            model = mat4::Identity();
        } else {
            model = bigModel;
        }
    }

    void draw(const mat4& view, GLuint pid){
        glBindVertexArray(_vao);
        ///--- Bind textures
        glActiveTexture(GL_TEXTURE4);
        glBindTexture(GL_TEXTURE_2D, _tex);
        glActiveTexture(GL_TEXTURE5);
        glBindTexture(GL_TEXTURE_2D, _texNormal);



        float texRatio = 6.0;

        glUniform1i(glGetUniformLocation(pid, "tex"), 4 /*GL_TEXTURE4*/);
        glUniform1i(glGetUniformLocation(pid, "normalTex"), 5 /*GL_TEXTURE5*/);
        glUniform1i(glGetUniformLocation(pid, "use_color"), GL_FALSE);
        // We don't use the normals for the ground or it becomes too dark
        glUniform1i(glGetUniformLocation(pid, "has_normal"), GL_FALSE);
        glUniform1f(glGetUniformLocation(pid, "texRatio"), texRatio);

        glUniformMatrix4fv(glGetUniformLocation(pid, "model"), 1, GL_FALSE, model.data());
        glUniformMatrix4fv(glGetUniformLocation(pid, "view"), 1, GL_FALSE, view.data());

        ///--- Draw
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        check_error_gl();
    }
};
