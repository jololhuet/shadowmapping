#pragma once
#include "icg_common.h"
#include <list>
#include "../attrib_locations.h"

bool rotate_cube = false;

class Mesh {
protected:
    GLuint _vao; ///< vertex array object
    vec3 _color; ///< Mesh Color
    opengp::Surface_mesh _mesh;

    Surface_mesh::Vertex_property<Point> vpoint;

public:
    void init(const char* obj = "tangle_cube.obj"){
        _color = vec3(1.0, 0.0, 0.0);

        ///--- Vertex one vertex Array
        glGenVertexArrays(1, &_vao);
        glBindVertexArray(_vao);

        _mesh.read(obj);
        _mesh.triangulate();
        _mesh.update_vertex_normals();

        ///--- Vertex coordinates
        {
            vpoint = _mesh.get_vertex_property<Point>("v:point");
            GLuint vbo;
            glGenBuffers(1, &vbo);
            glBindBuffer(GL_ARRAY_BUFFER, vbo);
            glBufferData(GL_ARRAY_BUFFER, _mesh.n_vertices() * sizeof(vec3), vpoint.data(), GL_STATIC_DRAW);

            ///--- Attribute
            glEnableVertexAttribArray(ATTRIB_LOC_vpoint);
            glVertexAttribPointer(ATTRIB_LOC_vpoint, 3, GL_FLOAT, DONT_NORMALIZE, ZERO_STRIDE, ZERO_BUFFER_OFFSET);
        }

        ///--- Normal coordinates
        {
            Surface_mesh::Vertex_property<Normal> vnormal = _mesh.get_vertex_property<Normal>("v:normal");
            GLuint vbo;
            glGenBuffers(1, &vbo);
            glBindBuffer(GL_ARRAY_BUFFER, vbo);
            glBufferData(GL_ARRAY_BUFFER, _mesh.n_vertices() * sizeof(vec3), vnormal.data(), GL_STATIC_DRAW);

            ///--- Attribute
            glEnableVertexAttribArray(ATTRIB_LOC_vnormal);
            glVertexAttribPointer(ATTRIB_LOC_vnormal, 3, GL_FLOAT, DONT_NORMALIZE, ZERO_STRIDE, ZERO_BUFFER_OFFSET);
        }


        ///--- Index Buffer
        {
            std::vector<unsigned int> indices;
            for(Surface_mesh::Face_iterator fit = _mesh.faces_begin(); fit != _mesh.faces_end(); ++fit) {
                unsigned int n = _mesh.valence(*fit);
                Surface_mesh::Vertex_around_face_circulator vit = _mesh.vertices(*fit);
                for(unsigned int v = 0; v < n; ++v) {
                    indices.push_back((*vit).idx());
                    ++vit;
                }
            }

            GLuint vbo;
            glGenBuffers(1, &vbo);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);
        }
    }

    void draw(const mat4& model, const mat4& view, GLuint pid){
        glBindVertexArray(_vao);

        glUniform3fv(glGetUniformLocation(pid, "mesh_color"), 1, _color.data());
        glUniform1i(glGetUniformLocation(pid, "use_color"), 1);
        glUniform1i(glGetUniformLocation(pid, "has_normal"), 0);

        glUniformMatrix4fv(glGetUniformLocation(pid, "model"), 1, GL_FALSE, model.data());
        glUniformMatrix4fv(glGetUniformLocation(pid, "view"), 1, GL_FALSE, view.data());

        ///--- Draw
        glDrawElements(GL_TRIANGLES, 3 * _mesh.n_faces(), GL_UNSIGNED_INT, ZERO_BUFFER_OFFSET);
    }

    void getVertices (std::list<std::vector<float> >& mesh_points, mat4& model) {
        Surface_mesh::Vertex_iterator vit;
        for (vit = _mesh.vertices_begin(); vit != _mesh.vertices_end(); ++vit)
        {
            vec3 vp = vpoint[*vit];
            vec4 vpoint_model = model * vec4(vp.x(), vp.y(), vp.z(), 1.0);

            std::vector<float> p(3);
            p[0] = vpoint_model.x();
            p[1] = vpoint_model.y();
            p[2] = vpoint_model.z();

            mesh_points.push_back(p);
        }
    }
};
