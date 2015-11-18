get_filename_component(DIR ${CMAKE_CURRENT_LIST_DIR} NAME)
file(GLOB_RECURSE SOURCES "*.cpp")
file(GLOB_RECURSE HEADERS "*.h")
file(GLOB_RECURSE SHADERS "*.glsl")
add_executable(${DIR} ${SOURCES} ${HEADERS} ${SHADERS})
target_link_libraries(${DIR} ${COMMON_LIBS})

# Deploying shaders
target_deploy_file(${DIR} _shadow_shaders/csm_vshader.glsl)
target_deploy_file(${DIR} _shadow_shaders/csm_fshader.glsl)
target_deploy_file(${DIR} _shadow_shaders/vsm_vshader.glsl)
target_deploy_file(${DIR} _shadow_shaders/vsm_fshader.glsl)
target_deploy_file(${DIR} _shadow_shaders/filter_vshader.glsl)
target_deploy_file(${DIR} _shadow_shaders/filter_fshader.glsl)
target_deploy_file(${DIR} _shadow_shaders/shadow_map_vshader.glsl)
target_deploy_file(${DIR} _shadow_shaders/shadow_map_fshader.glsl)
target_deploy_file(${DIR} _shadow_shaders/bias_vshader.glsl)
target_deploy_file(${DIR} _shadow_shaders/bias_fshader.glsl)
target_deploy_file(${DIR} _shadow_shaders/depth_vshader.glsl)
target_deploy_file(${DIR} _shadow_shaders/depth_fshader.glsl)


# Deploying textures
target_deploy_file(${DIR} _floor/planks.tga)
target_deploy_file(${DIR} _floor/planks_normal.tga)

# Deploying differents obj
target_deploy_file(${DIR} _mesh/tangle_cube.obj)
target_deploy_file(${DIR} _mesh/cube.obj)
target_deploy_file(${DIR} _mesh/sphere.obj)
target_deploy_file(${DIR} _mesh/tea.obj)
target_deploy_file(${DIR} _mesh/csm_cubes.obj)

target_deploy_file(${DIR} ../external/AntTweakBar/lib/libAntTweakBar.so)

