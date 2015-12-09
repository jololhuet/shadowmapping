#include "icg_common.h"
#include "ShadowBuffer.h"
#include "trackball.h"
#include "attrib_locations.h"
#include "_floor/Floor.h"
#include "_mesh/Mesh.h"
#include "Miniball.hpp"

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

// Fonctions:
void drawScene(bool _use_big_scene, GLuint _pid);

struct Frustum
{
    float near;
    float far;
    float aspect;
    float fovy;
    vec3 points[8];
    mat4 shadProj;
};

#define MAX_SPLITS 4
Frustum fs [MAX_SPLITS];
ShadowBuffer sb[MAX_SPLITS];  // FBO for shadow map generation

int width = 1024;
int height = 768;

// Scene elements
Mesh coco;
Mesh tea;
Mesh tangle_cube;
Mesh sphere;
Mesh wall;
Mesh csm_scene;
Floor ground_floor;

Mesh debug_cube;

GLuint default_pid;  // Handle for the default shader program
GLuint shadow_pid;  // Handle for the shadow map generation shader program

bool _use_big_scene = false;
bool _use_csm = false;

GLuint vsm_pid;     // Handle for the variance shadow mapping
GLuint bias_pid;    // Handle for the bias shadow map
GLuint filter_pid;  // Handle for the filtered shadow map
GLuint depth_pid;   // Handle to print the depths shadow map

bool use_polygon_offset = false;
bool use_normal_offset = false;
bool use_slope_bias = true;

GLuint depth_tex[MAX_SPLITS];  // Handle for the shadow map
GLuint depthVSM_tex[MAX_SPLITS];  // Handle for the shadow map

mat4 projection;  // Projection matrix for camera

mat4 view;  // View matrix for camera
mat4 view_big_scene;

mat4 trackball_matrix;
mat4 trackball_matrix_big_scene;

Trackball trackball;
int savedYMousePos = 0;

// Light stuff
bool _light_type = 0;
vec3 light_dir;  // Direction towards the light
vec3 light_pos;
mat4 light_projection;  // Projection matrix for light source
mat4 light_projection_orthogonal;
mat4 light_projection_perspective;
mat4 light_view;

mat4 offset_matrix;  // Affine transformation to map components from [-1, 1] to [0, 1], defined in init()

float bias = 0.005f;

int prev_buffer_size;
int buffer_size = 512;

float pcfSpread = 600.0;
float shadow_darkness = 0.20;

float vsmBlurRadius = 0.005;
int vsmBlurSteps = 4;

// Parameters for light perspective projection
GLfloat _light_near_pane = 1.0f;
GLfloat _light_far_pane = 50.0f;
GLfloat _light_fovy = 45.0f;
GLfloat _light_aspect = 1.0f;
GLfloat old_near = _light_near_pane;
GLfloat old_far = _light_far_pane;
GLfloat old_fovy = _light_fovy;
GLfloat old_aspect = _light_aspect;


//Miniball
typedef std::list<std::vector<float> >::const_iterator PointIterator;
typedef std::vector<float>::const_iterator CoordIterator;

typedef Miniball::Miniball <
        Miniball::CoordAccessor<PointIterator, CoordIterator> >
        MB;


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

    if (_use_big_scene) {
        view_big_scene = view_big_scene * trans;
    } else {
        view = view * trans;
    }

    keyboardState.lastUpdated = time;
}

vec3 lightCenterLargeScene;
vec3 lightCenterSmallScene;
float lightRadiusLargeScene;
float lightRadiusSmallScene;
mat4 createLightProjection (vec3 center, float radius, bool isOrtho){
    vec3 lightPosition;

    if (isOrtho) {
        lightPosition = light_dir;
    } else {
        lightPosition = light_pos;
    }

    float error = radius / 4.0;

    vec3 diff = center - lightPosition;
    float distance = sqrt(diff.dot(diff));

    float near = distance - radius -error;
    float far = distance + radius +error;

    if (isOrtho) {
        float mult = 1.4f;
//        cout << "near : " << near << endl;
//        cout << "far : " << far << endl;
        return OrthographicProjection(-mult*radius, mult*radius, -mult*radius, mult*radius, mult*near, mult*far);
    } else {
        float aspect = (float)width/(float)height;
        float fovy = 2 * asin(radius/distance) * 180.0f/ M_PI;
        _light_near_pane = near;
        _light_far_pane = far;
        _light_aspect = aspect;
        _light_fovy = fovy;

        return PerspectiveProjection(_light_fovy, _light_aspect, _light_near_pane, _light_far_pane);
    }
}


float camera_near = 0.1f;
float camera_far = 1000.0f;
void keyboard(int key, int action){
    vec3 center;
    float radius;
    if (_use_big_scene) {
        center = lightCenterLargeScene;
        radius = lightRadiusLargeScene;
    } else {
        center = lightCenterSmallScene;
        radius = lightRadiusSmallScene;
    }
     switch (key) {
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
            light_projection = createLightProjection(center, radius, true);
            _light_type = 0;
            bias = 0.005;
            std::cout<<"Directional Light"<<std::endl<<std::flush;
            break;

        case 'P':
            if(action != GLFW_RELEASE) return;
            light_projection = createLightProjection(center, radius, false);
            _light_type = 1;
            if (default_pid != vsm_pid) {
                bias = 0.10;
            } else {
                bias = 0.05;
            }
            std::cout<<"Spot Light"<<std::endl<<std::flush;
            break;

            /// 'I', 'J', 'K' and 'L' used to move the light source point.
            /// 'N' and 'M' to take the light up and down
        case 'L':
            if (_light_type == 0) light_dir += vec3(0.1f,0.0f,0.0f);
            else light_pos += vec3(0.3f,0.0f,0.0f);
            light_projection = createLightProjection(center, radius, _light_type == 0);
            break;

        case 'I':
            if (_light_type == 0) light_dir += vec3(0.0f,0.0f,-0.1f);
            else light_pos += vec3(0.0f,0.0f,-0.3f);
            light_projection = createLightProjection(center, radius, _light_type == 0);
            break;

        case 'J':
            if (_light_type == 0) light_dir += vec3(-0.1f,0.0f,0.0f);
            else light_pos += vec3(-0.3f,0.0f,0.0f);
            light_projection = createLightProjection(center, radius, _light_type == 0);
            break;

        case 'K':
            if (_light_type == 0) light_dir += vec3(0.0f,0.0f,0.1f);
            else light_pos += vec3(0.0f,0.0f,0.3f);
            light_projection = createLightProjection(center, radius, _light_type == 0);
            break;

        case 'N':
            if (_light_type == 0) light_dir += vec3(0.0f,0.1f,0.0f);
            else light_pos += vec3(0.0f,0.3f,0.0f);
            light_projection = createLightProjection(center, radius, _light_type == 0);
            break;

        case 'M':
            if (_light_type == 0) light_dir += vec3(0.0f,-0.1f,0.0f);
            else light_pos += vec3(0.0f,-0.3f,0.0f);
            light_projection = createLightProjection(center, radius, _light_type == 0);
            break;

        case '0':
            if(action != GLFW_RELEASE) return;
            default_pid = depth_pid;
//            _use_csm = false;
            std::cout<<"Mode DEPTH PRINT"<<std::endl<<std::flush;
            break;

        case '1':
            if(action != GLFW_RELEASE) return;
            default_pid = bias_pid;
            if (_light_type==1) {
                bias = 0.1;
            } else {
                bias = 0.005;
            }
            std::cout<<"Mode BIAS"<<std::endl<<std::flush;
            break;

        case '2':
            if(action != GLFW_RELEASE) return;
            default_pid = filter_pid;
            if (_light_type==1) {
                bias = 0.1;
            } else {
                bias = 0.005;
            }
            std::cout<<"Mode FILTER"<<std::endl<<std::flush;
            break;

        case '3':
            if(action != GLFW_RELEASE) return;
            default_pid = vsm_pid;
            if (_light_type==1) {
                bias = 0.05;
            } else {
                bias = 0.005;
            }

            std::cout<<"Mode VSM"<<std::endl<<std::flush;
            break;

        case '4':
            if(action != GLFW_RELEASE) return;
            _use_csm = !_use_csm;
            if (_use_csm) {
                std::cout<<"Mode CSM"<<std::endl<<std::flush;
            } else {
                std::cout<<"No more CSM"<<std::endl<<std::flush;
            }
            break;

        /// To switch the scene from default to big and reverse.
        case 'B':
            if(action != GLFW_RELEASE) return;
            _use_big_scene = !_use_big_scene;
            if (_use_big_scene) {
                std::cout<<"Large Scene"<<std::endl<<std::flush;
                light_projection = createLightProjection(lightCenterLargeScene, lightRadiusLargeScene, _light_type == 0);
                camera_near = 1.5*lightRadiusLargeScene;
                camera_far = 4.5*lightRadiusLargeScene;
                projection = PerspectiveProjection(45.0f, (GLfloat)width/height, camera_near, camera_far);
            } else {
                std::cout<<"Normal scene"<<std::endl<<std::flush;
                light_projection = createLightProjection(lightCenterSmallScene, lightRadiusSmallScene, _light_type == 0);
                camera_near = 1.5*lightRadiusSmallScene;
                camera_far = 4.5*lightRadiusSmallScene;
                projection = PerspectiveProjection(45.0f, (GLfloat)width/height, camera_near, camera_far);
            }
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

    projection = PerspectiveProjection(45.0f, (GLfloat)width / height, camera_near, camera_far);
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

// Compute the 8 corner points of the current view frustum
void updateFrustumPoints (Frustum& f, vec3 cam_pos, vec3 view_dir) {
    vec3 up = vec3(0.0, 1.0, 0.0);
    vec3 right = view_dir.cross(up);

    vec3 fc = cam_pos + view_dir*f.far;
    vec3 nc = cam_pos + view_dir*f.near;

    right = right.normalized();
    up = right.cross(view_dir);
    up = up.normalized();

    // these heights and widths are half the heights and widths of
    // the near and far plane rectangles
    float near_height = tan(f.fovy/2.0f) * f.near;
    float near_width = near_height * f.aspect;
    float far_height = tan(f.fovy/2.0f) * f.far;
    float far_width = far_height * f.aspect;

    f.points[0] = nc - up*near_height - right*near_width;
    f.points[1] = nc + up*near_height - right*near_width;
    f.points[2] = nc + up*near_height + right*near_width;
    f.points[3] = nc - up*near_height + right*near_width;

    f.points[4] = fc - up*far_height - right*far_width;
    f.points[5] = fc + up*far_height - right*far_width;
    f.points[6] = fc + up*far_height + right*far_width;
    f.points[7] = fc - up*far_height + right*far_width;

}

// Split the frustum into cur_num_splits smaller frustums
int cur_num_splits = 1;
void updateSplitDist ()  {
    float near = camera_near;
    float far = camera_far;

    float lambda = 0.75f;
    float ratio = far/near;
    fs[0].near = near;

    for(int i=1; i<cur_num_splits; i++)
    {
        float si = i / (float)cur_num_splits;

        fs[i].near = lambda*(near*pow(ratio, si)) + (1-lambda)*(near + (far - near)*si);
        fs[i-1].far = fs[i].near * 1.005f;
    }
    fs[cur_num_splits-1].far = far;
}

mat4 applyCropMatrix (Frustum& f, mat4 light_mv) {
    float minX = 1000;
    float maxX = -1000;
    float minY = 1000;
    float maxY = -1000;
    float minZ;
    float maxZ;

    vec4 transf = light_mv * vec4(f.points[0].x(), f.points[0].y(), f.points[0].z(), 1.0f);
    minZ = transf.z()/transf.w();
    maxZ = transf.z()/transf.w();

    for (int i = 0; i < 8; i++) {
        transf = light_mv * vec4(f.points[i].x(), f.points[i].y(), f.points[i].z(), 1.0f);
        if (transf.z() < minZ) minZ = transf.z();
        if (transf.z() > maxZ) maxZ = transf.z();
    }

    mat4 projection = OrthographicProjection(-1.0, 1.0, -1.0, 1.0, -maxZ, -minZ);

    mat4 light_mvp = projection * light_mv;

    for(int i=0; i<8; i++)
    {
        transf = light_mvp*vec4(f.points[i].x(), f.points[i].y(), f.points[i].z(), 1.0f);

        transf.x() /= transf.w();
        transf.y() /= transf.w();

        if(transf.x() > maxX) maxX = transf.x();
        if(transf.x() < minX) minX = transf.x();
        if(transf.y() > maxY) maxY = transf.y();
        if(transf.y() < minY) minY = transf.y();
    }

//    cout << "Min : " << minZ << " ; Max : " << maxZ << endl;

    float scaleX = 2.0f/(maxX - minX);
    float scaleY = 2.0f/(maxY - minY);
    float offsetX = -0.5f*(maxX + minX)*scaleX;
    float offsetY = -0.5f*(maxY + minY)*scaleY;

    light_mvp = mat4::Identity();
    light_mvp(0,0) = scaleX;
    light_mvp(1,1) = scaleY;
    light_mvp(0,3) = offsetX;
    light_mvp(1,3) = offsetY;
//    cout << light_mvp << endl;
//    cout << "-----------------------------------------" << endl;

    projection = light_mvp.transpose() * projection;

    return projection;
}

void drawDebugCube (GLfloat _pid, mat4& projection, mat4& view) {
    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    debug_cube.draw(projection, view, _pid);
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
}

mat4 all_mvp[MAX_SPLITS];
mat4 makeShadowMap (bool bigScene = false, bool _use_csm = false) {

    mat4 viewMat;
    vec3 center;
    float radius;

    mat4 light_mv = lookAt(vec3(0.0,0.0,0.0), vec3(-light_dir), vec3(-1.0f,0,0));
    mat4 depth_vp;

    if (!bigScene) {
        viewMat = view;
        center = lightCenterSmallScene;
        radius = lightRadiusSmallScene;
    } else {
        viewMat = view_big_scene;
        center = lightCenterLargeScene;
        radius = lightRadiusLargeScene;
    }

    vec3 cameraPosition = vec3(viewMat.col(3).x(), viewMat.col(3).y(), viewMat.col(3).z());
    vec3 view_dir = vec3(cameraPosition.normalized());

    // compute the z-distances for each split as seen in camera space
    updateSplitDist();

    if (_use_csm) {
        for (int i = 0; i < cur_num_splits; i++) {
            sb[i].bind();
                glActiveTexture(GL_TEXTURE0);
                glBindTexture(GL_TEXTURE_2D, depthVSM_tex[i]);
                check_error_gl();

                glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
                glViewport(0,0,buffer_size,buffer_size);

                // compute the camera frustum slice boundary points in world space
                updateFrustumPoints(fs[i], cameraPosition, view_dir);

                fs[i].shadProj = applyCropMatrix(fs[i], light_mv);
    //            cout << fs[i].shadProj << endl;
    //            cout << "-----------------------------------------" << endl;

                light_projection = fs[i].shadProj;
                light_view = light_mv;

                depth_vp = light_projection * light_view;
                glUniformMatrix4fv(glGetUniformLocation(shadow_pid, "depth_vp"), 1, GL_FALSE, depth_vp.data());

                // Draw the scene !
                drawScene(_use_csm, shadow_pid);
//                drawDebugCube(shadow_pid, light_projection, light_view);

                // We store the mvp matrices for each frustum for after
                all_mvp[i] = fs[i].shadProj * light_mv;
            sb[i].unbind();
        }
    } else {
        light_projection = createLightProjection(center, radius, _light_type == 0);
        sb[0].bind();
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, depthVSM_tex[0]);
            check_error_gl();

            glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
            glViewport(0,0,buffer_size,buffer_size);

            light_view = light_mv;
            depth_vp = light_projection * light_view;
            glUniformMatrix4fv(glGetUniformLocation(shadow_pid, "depth_vp"), 1, GL_FALSE, depth_vp.data());

            drawScene(_use_csm, shadow_pid);
//            drawDebugCube(shadow_pid, light_projection, light_view);
        sb[0].unbind();
    }

    return depth_vp;
}

mat4 createTransMatrix (float t1, float t2, float t3) {
    return Eigen::Affine3f(Eigen::Translation3f(vec3(t1, t2, t3))).matrix();
}

mat4 createScaleMatrix(float s1, float s2, float s3) {
    mat4 scale;
    scale <<    s1,  0.0f,  0.0f,  0.0f,
              0.0f,  s2,    0.0f,  0.0f,
              0.0f,  0.0f,  s3,    0.0f,
              0.0f,  0.0f,  0.0f,  1.0f;

    return scale;
}

mat4 generateModelMatrix (mat4& trans, mat4& scale) {
    return trans * scale;
}

mat4 cube_model_matrix;
mat4 sphere_model_matrix;
mat4 tea_model_matrix;
mat4 wall_model_matrix;

mat4 cone_model_matrix1;
mat4 cone_model_matrix2;
mat4 cone_model_matrix3;

mat4 csm_model_matrix;

void init() {
    cout << "Initialisation" << endl;
    // Light properties
    light_dir = vec3(-5.0, 10.0, 5.0);
    light_pos = vec3(light_dir);
    light_dir.normalize();


#ifdef WITH_ANTTWEAKBAR
    /*
     *  AntTweakBar stuff
     */
    {
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

        TwAddVarRW(bar, "Near Perspective", TW_TYPE_FLOAT, &_light_near_pane, " step=0.5 ");
        TwAddVarRW(bar, "Far Perspective", TW_TYPE_FLOAT, &_light_far_pane, " step=0.5 ");
        TwAddVarRW(bar, "Fovy Perspective", TW_TYPE_FLOAT, &_light_fovy, " step=0.5 ");
        TwAddVarRW(bar, "Aspect Perspective", TW_TYPE_FLOAT, &_light_aspect, " step=0.5 ");

        TwAddVarRW(bar, "VSM Blur Radius", TW_TYPE_FLOAT, &vsmBlurRadius, " step=0.0005 ");
        TwAddVarRW(bar, "VSM Blur #Steps", TW_TYPE_INT16, &vsmBlurSteps, " step=1 ");
    }
#endif

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);

    // Linking Shaders :
    {
        depth_pid = opengp::load_shaders("depth_vshader.glsl", "depth_fshader.glsl");
        glBindAttribLocation(depth_pid, ATTRIB_LOC_vpoint, "vpoint");
        glBindAttribLocation(depth_pid, ATTRIB_LOC_vtexcoord, "vtexcoord");
        glLinkProgram(depth_pid);

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
    }

    glViewport(0,0,width,height);

    ground_floor.init();

    // Default scene
    {
        // Meshes init
        coco.init("thin_cone.obj");
        tea.init("tea.obj");
        tangle_cube.init();
        sphere.init("sphere.obj");
        wall.init("cube.obj");
        debug_cube.init("big_cube.obj");

        mat4 cube_scale = createScaleMatrix(0.65, 0.65, 0.65);
        mat4 sphere_scale = createScaleMatrix(0.30, 0.30, 0.30);
        mat4 tea_scale = createScaleMatrix(0.10, 0.10, 0.10);
        mat4 wall_scale = createScaleMatrix(0.08, 0.80, 2.0);
        mat4 cone_scale = createScaleMatrix(0.2, 0.25, 0.2);

        mat4 cube_trans = createTransMatrix(0.0f, 0.33f, 0.0f);
        mat4 sphere_trans = createTransMatrix(-0.8f, 1.0f, 0.3f);
        mat4 tea_trans = createTransMatrix(-0.6f, 0.00f, -0.5f);
        mat4 wall_trans = createTransMatrix(-1.0f, 0.0f, -1.0f);

        mat4 cone_trans1 = createTransMatrix(0.1, 0.0, 0.8);
        mat4 cone_trans2 = createTransMatrix(0.4, 0.0, 0.6);
        mat4 cone_trans3 = createTransMatrix(0.7, 0.0, 0.8);

        cube_model_matrix = generateModelMatrix(cube_trans, cube_scale);
        sphere_model_matrix = generateModelMatrix(sphere_trans, sphere_scale);
        tea_model_matrix = generateModelMatrix(tea_trans, tea_scale);
        wall_model_matrix = generateModelMatrix(wall_trans, wall_scale);

        cone_model_matrix1 = generateModelMatrix(cone_trans1, cone_scale);
        cone_model_matrix2 = generateModelMatrix(cone_trans2, cone_scale);
        cone_model_matrix3 = generateModelMatrix(cone_trans3, cone_scale);

        std::list<std::vector<float> > all_points;
        tea.getVertices(all_points, tea_model_matrix);
        tangle_cube.getVertices(all_points, cube_model_matrix);
        sphere.getVertices(all_points, sphere_model_matrix);
        wall.getVertices(all_points, wall_model_matrix);
        coco.getVertices(all_points, cone_model_matrix1);
        coco.getVertices(all_points, cone_model_matrix2);
        coco.getVertices(all_points, cone_model_matrix3);

        MB mb (3, all_points.begin(), all_points.end());

        float ball_radius = sqrt(mb.squared_radius());
        vec3 center = vec3(mb.center());

        vec3 cameraPosition = center + vec3(0.0f, ball_radius, 3.0f*ball_radius);

        // Setting the projection matrix such that the whole scene is centered
        camera_near = 1.5*ball_radius;
        camera_far = 4.5*ball_radius;
        projection = PerspectiveProjection(45.0f, (GLfloat)width/height, camera_near, camera_far);

        view = lookAt(cameraPosition, vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));

        lightCenterSmallScene = center;
        lightRadiusSmallScene = ball_radius;
    }

    //Large Scene
    {
        // Mesh init
        csm_scene.init("csm_cubes.obj");

        // CSM scene
        float scsm = 1.0;
        mat4 csm_scale = createScaleMatrix(scsm, scsm, scsm);

        mat4 csm_trans = createTransMatrix(0.0f, 1.0f, 0.0f);

        csm_model_matrix = generateModelMatrix(csm_trans, csm_scale);

        std::list<std::vector<float> > all_points;
        csm_scene.getVertices(all_points, csm_model_matrix);

        MB mb (3, all_points.begin(), all_points.end());

        float ball_radius = sqrt(mb.squared_radius());
        vec3 center = vec3(mb.center());

        vec3 cameraPosition = center + vec3(0.0f, ball_radius, 3.0f*ball_radius);

        view_big_scene = lookAt(cameraPosition, vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));

        lightCenterLargeScene = center;
        lightRadiusLargeScene = ball_radius;
    }

    light_projection_orthogonal = createLightProjection(lightCenterSmallScene, lightRadiusSmallScene, true);
    light_projection_perspective = createLightProjection(lightCenterSmallScene, lightRadiusSmallScene, false);

    light_projection = light_projection_orthogonal;

    check_error_gl();

    default_pid = depth_pid;
    glUseProgram(default_pid);
    glUniform1i(glGetUniformLocation(default_pid, "shadow_map"), 1/*GL_TEXTURE1*/);

    default_pid = filter_pid;
    glUseProgram(default_pid);
    glUniform1i(glGetUniformLocation(default_pid, "shadow_map"), 1/*GL_TEXTURE1*/);

    default_pid = bias_pid;
    glUseProgram(default_pid);
    glUniform1i(glGetUniformLocation(default_pid, "shadow_map"), 1/*GL_TEXTURE1*/);

    default_pid = vsm_pid;
    glUseProgram(default_pid);
    glUniform1i(glGetUniformLocation(default_pid, "shadow_map"), 1/*GL_TEXTURE1*/);
    glUniform1i(glGetUniformLocation(default_pid, "vsm_shadow_map"), 0 /*GL_TEXTURE0*/);

    trackball_matrix = mat4::Identity();
    trackball_matrix_big_scene = mat4::Identity();

    // Matrix that can be used to move a point's components from [-1, 1] to [0, 1].
    offset_matrix << 0.5f, 0.0f, 0.0f, 0.5f,
                     0.0f, 0.5f, 0.0f, 0.5f,
                     0.0f, 0.0f, 0.5f, 0.5f,
                     0.0f, 0.0f, 0.0f, 1.0f;

    check_error_gl();

    for (int i = 0; i < cur_num_splits; i++) {
        sb[i].setSize(buffer_size, buffer_size);
        vec2 texs = sb[i].init();
        depth_tex[i] = texs[0];
        depthVSM_tex[i] = texs[1];
    }
    prev_buffer_size = buffer_size;

    check_error_gl();

}

// Draw the scene:
void drawScene (bool _use_big_scene, GLuint _pid) {
    if (!_use_big_scene) {
        wall.draw(wall_model_matrix, view * trackball_matrix, _pid);
        tangle_cube.draw(cube_model_matrix, view * trackball_matrix, _pid);
        sphere.draw(sphere_model_matrix, view * trackball_matrix, _pid);
        tea.draw(tea_model_matrix, view * trackball_matrix, _pid);
        ground_floor.setSize(0);
        ground_floor.draw(view * trackball_matrix, _pid);
        coco.draw(cone_model_matrix1, view * trackball_matrix, _pid);
        coco.draw(cone_model_matrix2, view * trackball_matrix, _pid);
        coco.draw(cone_model_matrix3, view * trackball_matrix, _pid);
    } else {
        csm_scene.draw(csm_model_matrix, view_big_scene*trackball_matrix_big_scene, _pid);
        ground_floor.setSize(1);
        ground_floor.draw(view_big_scene * trackball_matrix_big_scene, _pid);
    }
}

bool first = true;
void display() {
    if (first) {
//    first = false;
    check_error_gl();

    opengp::update_title_fps("Shadow Mapping");
    moveAll();

    glUseProgram(shadow_pid);
    check_error_gl();

    light_dir.normalize();

    // Check if we are using vsm shaders
    bool is_vsm = (default_pid == vsm_pid);

    if (buffer_size != prev_buffer_size) {
        prev_buffer_size = buffer_size;
        for (int i = 0; i < cur_num_splits; i++) {
            sb[i].setSize(buffer_size, buffer_size);
            vec2 texs = sb[i].init();
            depth_tex[i] = texs[0];
            depthVSM_tex[i] = texs[1];
        }
    }

    check_error_gl();

    if (_light_type == 0) {
        vec3 ortho_light_pos = light_dir;
        light_view = Eigen::lookAt(vec3(0,0,0), vec3(-ortho_light_pos), vec3(0,1,0));
    } else {
        light_view = Eigen::lookAt(light_pos,vec3(0,0,0),vec3(0,1,0));
        if (old_near != _light_near_pane || old_far != _light_far_pane || old_fovy != _light_fovy || old_aspect != _light_aspect) {
            light_projection = PerspectiveProjection(_light_fovy, _light_aspect, _light_near_pane, _light_far_pane);
            old_near = _light_near_pane;
            old_far = _light_far_pane;
            old_fovy = _light_fovy;
            old_aspect = _light_aspect;
        }
    }

    glUniform1f(glGetUniformLocation(shadow_pid, "bias"), (use_slope_bias)?bias:-1.0);
    glUniform1i(glGetUniformLocation(shadow_pid, "isVsm"), is_vsm);

    mat4 depth_vp = makeShadowMap(_use_big_scene, _use_csm);

    glUseProgram(default_pid);

    glUniformMatrix4fv(glGetUniformLocation(default_pid, "projection"), 1, GL_FALSE, projection.data());
    glUniform1i(glGetUniformLocation(default_pid, "light_type"), _light_type);
    if (_light_type == 0) {
        glUniform3fv(glGetUniformLocation(default_pid, "light_d"), 1, light_dir.data());
    } else {
        glUniform3fv(glGetUniformLocation(default_pid, "light_pos"), 1, light_pos.data());
    }

    check_error_gl();


    // Set matrix to transform from world space into NDC and then into [0, 1] ranges.
    mat4 depth_vp_offset = offset_matrix * depth_vp;

    glUniformMatrix4fv(glGetUniformLocation(default_pid, "depth_vp"), 1, GL_FALSE, depth_vp.data());
    glUniformMatrix4fv(glGetUniformLocation(default_pid, "depth_vp_offset"), 1, GL_FALSE, depth_vp_offset.data());

    glUniform1f(glGetUniformLocation(default_pid, "bias"), bias);
    glUniform1i(glGetUniformLocation(default_pid, "useCsm"), _use_csm);

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

    drawScene(_use_big_scene, default_pid);

    glIsQuery(0);

    check_error_gl();

#ifdef WITH_ANTTWEAKBAR
    TwDraw();
#endif
    }
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
        if (_use_big_scene){
            old_trackball_matrix_CSM = trackball_matrix_big_scene;  // Store the current state of the model matrix.
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

        if (_use_big_scene){
            trackball_matrix_big_scene = trackball.drag(p.x(), p.y()) * old_trackball_matrix_CSM;
        } else {
            trackball_matrix = trackball.drag(p.x(), p.y()) * old_trackball_matrix;
        }
    }

    // Zoom
    if (glfwGetMouseButton(GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
        mat4 zRot = mat4::Identity();
        zRot(2, 3) = (y - savedYMousePos) / 50.0f;

        if (_use_big_scene) {
            view_big_scene = zRot * view_big_scene;
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
