#include "icg_common.h"
#include "ShadowBuffer.h"
#include "trackball.h"
#include "attrib_locations.h"
#include "_floor/Floor.h"
#include "_mesh/Mesh.h"


/**
 *  Precision, the code implemented here come from multiple differents tutorial and courses.
 *  The Introduction to computer graphics course code is re-used  and for the differents
 *  shadow mapping techniques we get inspired by multiples tutos, links to thoses are here :
 *
 *  http://www.opengl-tutorial.org/intermediate-tutorials/tutorial-16-shadow-mapping/
 *  http://www.dissidentlogic.com/old/
 *
 */
#ifdef WITH_ANTTWEAKBAR
#include "../external/AntTweakBar/include/AntTweakBar.h"
TwBar *bar;
#endif

int width = 1024;
int height = 768;

ShadowBuffer sb;  // FBO for shadow map generation

// Shadow Buffers for the cascaded shadow map
ShadowBuffer sb_near;
ShadowBuffer sb_far;

// Scene elements
Mesh tea;
Mesh tangle_cube;
Mesh sphere;
Floor ground_floor;
Mesh wall;
Mesh csm_scene;

//Light light;


GLuint default_pid;  // Handle for the default shader program
GLuint shadow_pid;  // Handle for the shadow map generation shader program

GLuint csm_pid;     // Handle for the cascaded shadow mapping
GLuint vsm_pid;     // Handle for the variance shadow mapping
GLuint bias_pid;    // Handle for the bias shadow map
GLuint filter_pid;  // Handle for the filtered shadow map
GLuint depth_pid;   // Handle to print the depths shadow map

bool use_polygon_offset = false;
bool use_normal_offset = false;
bool use_slope_bias = true;

bool isCSM = false;

GLuint depth_tex;  // Handle for the shadow map
GLuint depthVSM_tex;  // Handle for the shadow map

mat4 projection;  // Projection matrix for camera
mat4 view;  // View matrix for camera
mat4 view_CSM;  // View matrix for camera

mat4 trackball_matrix;
mat4 trackball_matrix_CSM;

Trackball trackball;
int savedYMousePos = 0;

bool _light_type = 0;
vec3 light_dir;  // Direction towards the light
vec3 light_pos;
mat4 light_projection;  // Projection matrix for light source

mat4 offset_matrix;  // Affine transformation to map components from [-1, 1] to [0, 1], defined in init()

float bias = 0.01f;
int prev_buffer_size;
int buffer_size = 1024;

float pcfSpread = 700.0;
float shadow_darkness = 0.15;

float vsmBlurRadius = 0.01;
int vsmBlurSteps = 4;

// Parameters for perspective projection
GLfloat _near_pane = 1.0f;
GLfloat _far_pane = 50.0f;
GLfloat _fovy = 45.0f;
GLfloat _aspect = 1.0f;
GLfloat old_near = _near_pane;
GLfloat old_far = _far_pane;
GLfloat old_fovy = _fovy;
GLfloat old_aspect = _aspect;


mat4 PerspectiveProjection(float fovy, float aspect, float near, float far){
    float fovyDegrees = fovy  * (M_PI / 180.0f);
    float halfHeight =  near * tanf(fovyDegrees / 2.0f);
    float top = halfHeight;
    float bottom = -top;
    float right = aspect * halfHeight;
    float left = -right;

    mat4 projection = mat4::Zero();
    projection(0, 0) = 2.0f * near / (right - left);
    projection(1, 1) = 2.0f * near / (top - bottom);
    projection(2, 2) = -(far + near) / (far - near);
    projection(0, 2) = (right + left) / (right - left);
    projection(1, 2) = (top + bottom) / (top - bottom);
    projection(2, 3) = -(2.0f * far * near) / (far - near);
    projection(3, 2) = -1;

    return projection;
}

mat4 OrthographicProjection(float left, float right, float bottom, float top, float near, float far){
    assert(right > left);
    assert(far > near);
    assert(top > bottom);
    mat4 ortho = mat4::Zero();
    ortho(0, 0) = 2.0f / (right - left);
    ortho(1, 1) = 2.0f / (top - bottom);
    ortho(2, 2) = -2.0f / (far - near);
    ortho(3, 3) = 1.0f;
    ortho(0, 3) = -(right + left) / (right - left);
    ortho(1, 3) = -(top + bottom) / (top - bottom);
    ortho(2, 3) = -(far + near) / (far - near);
    return ortho;
}

struct keyInput {
    // Keys to move the camera
    bool isWPressed;
    bool isAPressed;
    bool isDPressed;
    bool isSPressed;

    //Keys to move the light source
    bool isUpPressed;
    bool isDownPressed;
    bool isLeftPressed;
    bool isRightPressed;


    float lastUpdated;
    bool needUpdate;
};
struct keyInput keyboardState;


void moveAll() {
    const float velocity = 1.0f;
    float time = glfwGetTime();
    float d_time = time - keyboardState.lastUpdated;

    mat4 trans = mat4::Identity();

    if (keyboardState.isAPressed) {
        trans = trans * Eigen::Affine3f(Eigen::Translation3f(vec3(velocity*d_time, 0.0f, 0.0f))).matrix();
    }
    if (keyboardState.isSPressed) {
        trans = trans * Eigen::Affine3f(Eigen::Translation3f(vec3(0.0f, 0.0f, -velocity*d_time))).matrix();
    }
    if (keyboardState.isDPressed) {
        trans = trans * Eigen::Affine3f(Eigen::Translation3f(vec3(-velocity*d_time, 0.0f, 0.0f))).matrix();
    }
    if (keyboardState.isWPressed) {
        trans = trans * Eigen::Affine3f(Eigen::Translation3f(vec3(0.0f, 0.0f, velocity*d_time))).matrix();
    }

    if (isCSM) {
        view_CSM = view_CSM * trans;
    } else {
        view = view * trans;
    }

    keyboardState.lastUpdated = time;

}


void keyboard(int key, int action){
     switch (key) {
//        if(action != GLFW_RELEASE) return;

            /// 'W' 'A' 'S' and 'D' are use to move the camera around the scene
        case 'W':
            keyboardState.isWPressed = !keyboardState.isWPressed;
            break;

        case 'A':
            keyboardState.isAPressed = !keyboardState.isAPressed;
            break;

        case 'S':
            keyboardState.isSPressed = !keyboardState.isSPressed;
            break;

        case 'D':
            keyboardState.isDPressed = !keyboardState.isDPressed;
            break;

            /// 'P' and 'O' are used to change the type of light
        case 'O':
            if(action != GLFW_RELEASE) return;
            light_projection = OrthographicProjection(-1.5,1.5,-1.5,1.5,-0.5,2.5);
            _light_type = 0;
            std::cout<<"Directional Light"<<std::endl<<std::flush;
            break;

        case 'P':
            if(action != GLFW_RELEASE) return;
            light_projection = PerspectiveProjection(_fovy, _aspect, _near_pane, _far_pane);
            _light_type = 1;
            std::cout<<"Spot Light"<<std::endl<<std::flush;
            break;

            /// 'I', 'J', 'K' and 'L' used to move the light source point.
            /// 'N' and 'M' to take the light up and down
        case 'L':
            if (_light_type == 0) light_dir += vec3(0.1f,0.0f,0.0f);
            else light_pos += vec3(0.3f,0.0f,0.0f);
            break;

        case 'I':
            if (_light_type == 0) light_dir += vec3(0.0f,0.0f,-0.1f);
            else light_pos += vec3(0.0f,0.0f,-0.3f);
            break;

        case 'J':
            if (_light_type == 0) light_dir += vec3(-0.1f,0.0f,0.0f);
            else light_pos += vec3(-0.3f,0.0f,0.0f);
            break;

        case 'K':
            if (_light_type == 0) light_dir += vec3(0.0f,0.0f,0.1f);
            else light_pos += vec3(0.0f,0.0f,0.3f);
            break;

        case 'N':
            if (_light_type == 0) light_dir += vec3(0.0f,0.1f,0.0f);
            else light_pos += vec3(0.0f,0.3f,0.0f);
            break;

        case 'M':
            if (_light_type == 0) light_dir += vec3(0.0f,-0.1f,0.0f);
            else light_pos += vec3(0.0f,-0.3f,0.0f);
            break;

        case '0':
            if(action != GLFW_RELEASE) return;
            default_pid = depth_pid;
            std::cout<<"Mode DEPTH PRINT"<<std::endl<<std::flush;
            break;

        case '1':
            if(action != GLFW_RELEASE) return;
            default_pid = bias_pid;
            std::cout<<"Mode BIAS"<<std::endl<<std::flush;
            break;

        case '2':
            if(action != GLFW_RELEASE) return;
            default_pid = filter_pid;
            std::cout<<"Mode FILTER"<<std::endl<<std::flush;
            break;

        case '3':
            if(action != GLFW_RELEASE) return;
            default_pid = vsm_pid;
            std::cout<<"Mode VSM"<<std::endl<<std::flush;
            break;

        case '4':
            if(action != GLFW_RELEASE) return;
            default_pid = csm_pid;
            std::cout<<"Mode CSM"<<std::endl<<std::flush;
            break;

        default:
            break;
     }
}

// Gets called when the windows is resized.
void resize_callback(int w, int h) {
    width = w;
    height = h;

    std::cout << "Window has been resized to " << width << "x" << height << "." << std::endl;
    glViewport(0, 0, width, height);

    projection = PerspectiveProjection(45.0f, (GLfloat)width / height, 0.1f, 100.0f);
}

void TW_CALL SlopeBiasSet (const void *value, void *clientData)
{
    use_slope_bias = *(const bool *)value;  // for instance
    if (use_polygon_offset) {
        use_polygon_offset = false;
    }
    if (use_normal_offset) {
        use_normal_offset = false;
    }
}
void TW_CALL SlopeBiasGet (void *value, void *clientData)
{
    *(bool *)value = use_slope_bias;  // for instance
}
void TW_CALL PolygonOffsetSet (const void *value, void *clientData)
{
    use_polygon_offset = *(const bool *)value;  // for instance
    if (use_slope_bias) {
        use_slope_bias = false;
    }
    if (use_normal_offset) {
        use_normal_offset = false;
    }
}
void TW_CALL PolygonOffsetGet (void *value, void *clientData)
{
    *(bool *)value = use_polygon_offset;  // for instance
}
void TW_CALL NormalOffsetSet (const void *value, void *clientData)
{
    use_normal_offset = *(const bool *)value;  // for instance
    if (use_polygon_offset) {
        use_polygon_offset = false;
    }
    if (use_slope_bias) {
        use_slope_bias = false;
    }
}
void TW_CALL NormalOffsetGet (void *value, void *clientData)
{
    *(bool *)value = use_normal_offset;  // for instance
}


void init() {
    light_dir = vec3(-5.0, 10.0, 5.0);
    light_pos = vec3(light_dir); //vec3(2.0, 3.0, 0.0);
    light_dir.normalize();

#ifdef WITH_ANTTWEAKBAR
    TwInit(TW_OPENGL_CORE, NULL);
    TwWindowSize(width, height);
    bar = TwNewBar("AntTweakBar");
    TwAddVarRW(bar, "LightDir", TW_TYPE_DIR3F, &light_dir, NULL);
    TwAddVarRW(bar, "Bias", TW_TYPE_FLOAT, &bias, " step=0.005 ");
    TwAddVarRW(bar, "Buffer size", TW_TYPE_INT32, &buffer_size, " step=512 ");

    TwAddVarCB(bar, "Use Slope Bias (default)", TW_TYPE_BOOLCPP, SlopeBiasSet, SlopeBiasGet, NULL, NULL);
    TwAddVarCB(bar, "Use Polygon Offset", TW_TYPE_BOOLCPP, PolygonOffsetSet, PolygonOffsetGet, NULL, NULL);
    TwAddVarCB(bar, "Use Normal Offset", TW_TYPE_BOOLCPP, NormalOffsetSet, NormalOffsetGet, NULL, NULL);

    TwAddVarRW(bar, "Pcf Spread", TW_TYPE_FLOAT, &pcfSpread, " step=50.0 ");
    TwAddVarRW(bar, "Shadow Darkness", TW_TYPE_FLOAT, &shadow_darkness, " step=0.05 ");

    TwAddVarRW(bar, "Near Perspective", TW_TYPE_FLOAT, &_near_pane, " step=0.5 ");
    TwAddVarRW(bar, "Far Perspective", TW_TYPE_FLOAT, &_far_pane, " step=0.5 ");
    TwAddVarRW(bar, "Fovy Perspective", TW_TYPE_FLOAT, &_fovy, " step=0.5 ");
    TwAddVarRW(bar, "Aspect Perspective", TW_TYPE_FLOAT, &_aspect, " step=0.5 ");

    TwAddVarRW(bar, "VSM Blur Radius", TW_TYPE_FLOAT, &vsmBlurRadius, " step=0.0005 ");
    TwAddVarRW(bar, "VSM Blur #Steps", TW_TYPE_INT16, &vsmBlurSteps, " step=1 ");

#endif

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);

    // Linking Shaders :
    depth_pid = opengp::load_shaders("depth_vshader.glsl", "depth_fshader.glsl");
    glBindAttribLocation(depth_pid, ATTRIB_LOC_vpoint, "vpoint");
    glBindAttribLocation(depth_pid, ATTRIB_LOC_vtexcoord, "vtexcoord");
    glLinkProgram(depth_pid);

    csm_pid = opengp::load_shaders("csm_vshader.glsl", "csm_fshader.glsl");
    glBindAttribLocation(csm_pid, ATTRIB_LOC_vpoint, "vpoint");
    glBindAttribLocation(csm_pid, ATTRIB_LOC_vtexcoord, "vtexcoord");
    glLinkProgram(csm_pid);

    vsm_pid = opengp::load_shaders("vsm_vshader.glsl", "vsm_fshader.glsl");
    glBindAttribLocation(vsm_pid, ATTRIB_LOC_vpoint, "vpoint");
    glBindAttribLocation(vsm_pid, ATTRIB_LOC_vtexcoord, "vtexcoord");
    glLinkProgram(vsm_pid);

    filter_pid = opengp::load_shaders("filter_vshader.glsl", "filter_fshader.glsl");
    glBindAttribLocation(filter_pid, ATTRIB_LOC_vpoint, "vpoint");
    glBindAttribLocation(filter_pid, ATTRIB_LOC_vtexcoord, "vtexcoord");
    glLinkProgram(filter_pid);

    bias_pid = opengp::load_shaders("bias_vshader.glsl", "bias_fshader.glsl");
    glBindAttribLocation(bias_pid, ATTRIB_LOC_vpoint, "vpoint");
    glBindAttribLocation(bias_pid, ATTRIB_LOC_vtexcoord, "vtexcoord");
    glLinkProgram(bias_pid);


    shadow_pid = opengp::load_shaders("shadow_map_vshader.glsl", "shadow_map_fshader.glsl");
    glBindAttribLocation(shadow_pid, ATTRIB_LOC_vpoint, "vpoint");
    glLinkProgram(shadow_pid);

    glViewport(0,0,width,height);

    float ratio = width / (float) height;
    projection = PerspectiveProjection(45.0f, ratio, 0.1f, 100.0f);
    view = Eigen::Affine3f(Eigen::Translation3f(0.0f, 0.1f, -4.0f)).matrix();
    view_CSM = Eigen::Affine3f(Eigen::Translation3f(0.0f, -0.5f, -5.0f)).matrix();

    default_pid = csm_pid;
    glUseProgram(default_pid);
    glUniformMatrix4fv(glGetUniformLocation(default_pid, "projection"), 1, GL_FALSE, projection.data());
    glUniform1i(glGetUniformLocation(default_pid, "shadow_map"), 1/*GL_TEXTURE1*/);

    default_pid = depth_pid;
    glUseProgram(default_pid);
    glUniformMatrix4fv(glGetUniformLocation(default_pid, "projection"), 1, GL_FALSE, projection.data());
    glUniform1i(glGetUniformLocation(default_pid, "shadow_map"), 1/*GL_TEXTURE1*/);

    default_pid = filter_pid;
    glUseProgram(default_pid);
    glUniformMatrix4fv(glGetUniformLocation(default_pid, "projection"), 1, GL_FALSE, projection.data());
    glUniform1i(glGetUniformLocation(default_pid, "shadow_map"), 1/*GL_TEXTURE1*/);

    default_pid = bias_pid;
    glUseProgram(default_pid);
    glUniformMatrix4fv(glGetUniformLocation(default_pid, "projection"), 1, GL_FALSE, projection.data());
    glUniform1i(glGetUniformLocation(default_pid, "shadow_map"), 1/*GL_TEXTURE1*/);

    default_pid = vsm_pid;
    glUseProgram(default_pid);
    glUniformMatrix4fv(glGetUniformLocation(default_pid, "projection"), 1, GL_FALSE, projection.data());
    glUniform1i(glGetUniformLocation(default_pid, "shadow_map"), 1/*GL_TEXTURE1*/);
    glUniform1i(glGetUniformLocation(default_pid, "vsm_shadow_map"), 0 /*GL_TEXTURE0*/);

    light_projection = OrthographicProjection(-1.5,1.5,-1.5,1.5,-0.5,2.5);

    trackball_matrix = mat4::Identity();
    trackball_matrix_CSM = mat4::Identity();

    // Matrix that can be used to move a point's components from [-1, 1] to [0, 1].
    offset_matrix << 0.5f, 0.0f, 0.0f, 0.5f,
                     0.0f, 0.5f, 0.0f, 0.5f,
                     0.0f, 0.0f, 0.5f, 0.5f,
                     0.0f, 0.0f, 0.0f, 1.0f;

    check_error_gl();


    csm_scene.init("csm_cubes.obj");
    tea.init("tea.obj");
    tangle_cube.init();
    sphere.init("sphere.obj");

    sb.setSize(buffer_size, buffer_size);
    prev_buffer_size = buffer_size;
    vec2 texs = sb.init();
    depth_tex = texs[0];
    depthVSM_tex = texs[1];

    check_error_gl();

    ground_floor.init();
    wall.init("cube.obj");
}

void makeTightFrustum () {
    //Do not know yet -- Preprocessing
}

// Split the frustum into 2 smaller frustums
void splitFrustum ()  {
    float near = 0;
    float far = 1.0;

    vec2 splitFar = vec2(0,0);
    vec2 splitNear = vec2(0,0);

    float lambda = 0.8;

    for (int i = 1; i <= 2; i++) {
        splitFar[i-1] = lambda * near*pow(far/near, i/2.0) + (1-lambda) * near + (far-near)*i/2.0;
    }
}

void display() {
    check_error_gl();

    opengp::update_title_fps("Shadow Mapping");
    moveAll();

    glUseProgram(shadow_pid);
    check_error_gl();

    // Default scene
    mat4 cube_scale;
    mat4 sphere_scale;
    mat4 tea_scale;
    mat4 wall_scale;

    mat4 cube_trans;
    mat4 sphere_trans;
    mat4 tea_trans;
    mat4 wall_trans;

    mat4 cube_model_matrix;
    mat4 sphere_model_matrix;
    mat4 tea_model_matrix;
    mat4 wall_model_matrix;

    // CSM scene
    mat4 csm_scale;

    mat4 csm_trans;

    mat4 csm_model_matrix;

    isCSM = default_pid == csm_pid;
    if (!isCSM) {
    // Scaling matrix to scale the differents mesh down to a reasonable size.
    float sc = 0.65;
    cube_scale << sc,    0.0f,  0.0f,  0.0f,
                  0.0f,  sc,    0.0f,  0.0f,
                  0.0f,  0.0f,  sc,    0.0f,
                  0.0f,  0.0f,  0.0f,  1.0f;

    float ss = 0.30;
    sphere_scale << ss,    0.0f,  0.0f,  0.0f,
                    0.0f,  ss,    0.0f,  0.0f,
                    0.0f,  0.0f,  ss,    0.0f,
                    0.0f,  0.0f,  0.0f,  1.0f;

    float st = 0.1;
    tea_scale <<    st,    0.0f,  0.0f,  0.0f,
                    0.0f,  st,    0.0f,  0.0f,
                    0.0f,  0.0f,  st,    0.0f,
                    0.0f,  0.0f,  0.0f,  1.0f;

    float sw = .8;
    wall_scale << 0.08,    0.0f,  0.0f,  0.0f,
                  0.0f,  sw,    0.0f,  0.0f,
                  0.0f,  0.0f,  2,    0.0f,
                  0.0f,  0.0f,  0.0f,  1.0f;




    // Transformations for the mesh models_matrixs
    cube_trans = Eigen::Affine3f(Eigen::Translation3f(vec3(0.0f, 0.33f, 0.0f))).matrix();
    sphere_trans = Eigen::Affine3f(Eigen::Translation3f(vec3(-0.8f, 1.0f, 0.3f))).matrix();
    tea_trans = Eigen::Affine3f(Eigen::Translation3f(vec3(-0.6f, 0.00f, -0.5f))).matrix();
    wall_trans = Eigen::Affine3f(Eigen::Translation3f(vec3(-1.0f, 0.0f, -1.0f))).matrix();

    cube_model_matrix = cube_trans * cube_scale;
    sphere_model_matrix = sphere_trans * sphere_scale;
    tea_model_matrix = tea_trans * tea_scale;
    wall_model_matrix = wall_trans * wall_scale;

    } else {
    float scsm = 0.5;
    csm_scale << scsm,    0.0f,  0.0f,  0.0f,
                  0.0f,  scsm,    0.0f,  0.0f,
                  0.0f,  0.0f,  scsm,    0.0f,
                  0.0f,  0.0f,  0.0f,  1.0f;

    csm_trans = Eigen::Affine3f(Eigen::Translation3f(vec3(0.0f, 0.5f, 0.0f))).matrix();

    csm_model_matrix = csm_trans * csm_scale;
    }

    light_dir.normalize();

    // Check if we are using vsm shaders
    bool is_vsm = (default_pid == vsm_pid);

    if (buffer_size != prev_buffer_size || is_vsm != sb.getIsVsm()) {
        sb.setVsm(is_vsm);
        sb.setSize(buffer_size, buffer_size);
        vec2 texs = sb.init();
        depth_tex = texs[0];
        depthVSM_tex = texs[1];
        prev_buffer_size = buffer_size;
    }

    check_error_gl();

//    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

    //=== BIND
    sb.bind();
    mat4 light_view;
    if (_light_type == 0) {
        light_view = Eigen::lookAt(light_dir,vec3(0,0,0),vec3(0,1,0));
    } else {
        light_view = Eigen::lookAt(light_pos,vec3(0,0,0),vec3(0,1,0));
        if (old_near != _near_pane || old_far != _far_pane || old_fovy != _fovy || old_aspect != _aspect) {
            light_projection = PerspectiveProjection(_fovy, _aspect, _near_pane, _far_pane);
            old_near = _near_pane;
            old_far = _far_pane;
            old_fovy = _fovy;
            old_aspect = _aspect;
        }
    }

    mat4 depth_vp = light_projection * light_view;
    glUniformMatrix4fv(glGetUniformLocation(shadow_pid, "depth_vp"), 1, GL_FALSE, depth_vp.data());
    glUniform1i(glGetUniformLocation(shadow_pid, "isVsm"), is_vsm);

    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    glCullFace(GL_BACK);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, depthVSM_tex);
    check_error_gl();




    if (!isCSM) {
        wall.draw(wall_model_matrix, view * trackball_matrix, shadow_pid);
        tangle_cube.draw(cube_model_matrix, view * trackball_matrix, shadow_pid);
        sphere.draw(sphere_model_matrix, view * trackball_matrix, shadow_pid);
        tea.draw(tea_model_matrix, view * trackball_matrix, shadow_pid);
        ground_floor.draw(view * trackball_matrix, shadow_pid);
    } else {
        csm_scene.draw(csm_model_matrix, view_CSM*trackball_matrix_CSM, shadow_pid);
        ground_floor.draw(view_CSM * trackball_matrix_CSM, shadow_pid);
    }


    sb.unbind();
    //=== UNBIND

    glUseProgram(default_pid);
    glUniform1i(glGetUniformLocation(default_pid, "light_type"), _light_type);
    if (_light_type == 0) {
        glUniform3fv(glGetUniformLocation(default_pid, "light_d"), 1, light_dir.data());
    } else {
        glUniform3fv(glGetUniformLocation(default_pid, "light_pos"), 1, light_pos.data());
    }

    check_error_gl();


    // Set matrix to transform from world space into NDC and then into [0, 1] ranges.
    mat4 depth_vp_offset = offset_matrix * depth_vp;


    glUniformMatrix4fv(glGetUniformLocation(default_pid, "light_projection"), 1, GL_FALSE, light_projection.data());
    glUniformMatrix4fv(glGetUniformLocation(default_pid, "light_view"), 1, GL_FALSE, light_view.data());
    glUniformMatrix4fv(glGetUniformLocation(default_pid, "depth_vp"), 1, GL_FALSE, depth_vp.data());
    glUniformMatrix4fv(glGetUniformLocation(default_pid, "depth_vp_offset"), 1, GL_FALSE, depth_vp_offset.data());

    glUniform1f(glGetUniformLocation(default_pid, "bias"), bias);

    glUniform1i(glGetUniformLocation(default_pid, "usePolygonOffset"), use_polygon_offset);
    glUniform1i(glGetUniformLocation(default_pid, "useNormalOffset"), use_normal_offset);
    glUniform1i(glGetUniformLocation(default_pid, "useSlopeBias"), use_slope_bias);

    glUniform1f(glGetUniformLocation(default_pid, "pcfSpread"), pcfSpread);
    glUniform1f(glGetUniformLocation(default_pid, "shadowDarkness"), shadow_darkness);

    // VSM uniforms :
    if (is_vsm) {
        glUniform1f(glGetUniformLocation(default_pid, "blurRadius"), vsmBlurRadius);
        glUniform1i(glGetUniformLocation(default_pid, "blurSteps"), vsmBlurSteps);
    }

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glCullFace(GL_FRONT);

    if (use_polygon_offset && !is_vsm) {
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(5, 0);
    } else {
        glDisable(GL_POLYGON_OFFSET_FILL);
    }




    if (!isCSM) {
        wall.draw(wall_model_matrix, view * trackball_matrix, default_pid);
        tangle_cube.draw(cube_model_matrix, view * trackball_matrix, default_pid);
        sphere.draw(sphere_model_matrix, view * trackball_matrix, default_pid);
        tea.draw(tea_model_matrix, view * trackball_matrix, default_pid);
        ground_floor.setSize(0);
        ground_floor.draw(view * trackball_matrix, default_pid);
    } else {
        csm_scene.draw(csm_model_matrix, view_CSM*trackball_matrix_CSM, default_pid);
        ground_floor.setSize(1);
        ground_floor.draw(view_CSM * trackball_matrix_CSM, default_pid);
    }


    glIsQuery(0);

    check_error_gl();

#ifdef WITH_ANTTWEAKBAR
    TwDraw();
#endif
}

/// Trackball things
// Transforms glfw screen coordinates into normalized OpenGL coordinates.
vec2 transform_screen_coords(int x, int y) {
    return vec2(2.0f * (float)x / width - 1.0f,
                1.0f - 2.0f * (float)y / height);
}

mat4 old_trackball_matrix;
mat4 old_trackball_matrix_CSM;

void mouse_button(int button, int action) {
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
        int x_i, y_i;
        glfwGetMousePos(&x_i, &y_i);
        vec2 p = transform_screen_coords(x_i, y_i);
        trackball.begin_drag(p.x(), p.y());
        if (isCSM){
            old_trackball_matrix_CSM = trackball_matrix_CSM;  // Store the current state of the model matrix.
        } else {
            old_trackball_matrix = trackball_matrix;  // Store the current state of the model matrix.
        }
    }

    if(button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
        glfwGetMousePos(NULL, &savedYMousePos);
    }
}

void mouse_pos(int x, int y) {
    if (glfwGetMouseButton(GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
        vec2 p = transform_screen_coords(x, y);

        if (isCSM){
            trackball_matrix_CSM = trackball.drag(p.x(), p.y()) * old_trackball_matrix_CSM;
        } else {
            trackball_matrix = trackball.drag(p.x(), p.y()) * old_trackball_matrix;
        }
    }

    // Zoom
    if (glfwGetMouseButton(GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
        mat4 zRot = mat4::Identity();
        zRot(2, 3) = (y - savedYMousePos) / 50.0f;

        if (isCSM) {
            view_CSM = zRot * view_CSM;
        } else {
            view = zRot * view;
        }

        glfwGetMousePos(NULL, &savedYMousePos);
    }
}

/// End of trackball things

// Mouse capture. By AntTweakBar for some things and by the program itself otherwise
void GLFWCALL OnMousePos (int mouseX, int mouseY) {
    if( !TwEventMousePosGLFW(mouseX, mouseY) ) {
        mouse_pos(mouseX, mouseY);
    }
}

// Keyboard capture. By AntTweakBar for some things and by the program itself otherwise
void GLFWCALL OnKeyClic (int button, int action) {
    if( !TwEventCharGLFW(button, action) && !TwEventKeyGLFW(button, action)) {
        keyboard(button, action);
    }
}

// Mouse capture. By AntTweakBar for some things and by the program itself otherwise
void GLFWCALL OnMouseButton (int button, int action) {
    if( !TwEventMouseButtonGLFW(button, action) ) {
        mouse_button(button, action);
    }
}


int main(int, char**) {
    glfwInitWindowSize(width, height);
    glfwCreateWindow("Shadow Mapping");
    glfwDisplayFunc(display);
    glfwSetWindowSizeCallback(&resize_callback);
    glfwSetKeyCallback(OnKeyClic);
    glfwSetMouseButtonCallback(OnMouseButton);
    glfwSetMousePosCallback(OnMousePos);
    init();
    keyboard(GLFW_KEY_KP_1, 0);
    glfwMainLoop();

#ifdef WITH_ANTTWEAKBAR
    TwTerminate();
#endif

    return EXIT_SUCCESS;    
}
