#version 330

in vec4 o_vertex;
in vec3 o_normal;
in vec2 o_texcoord;
in float o_material;

out vec3 v_vertex;
out vec3 v_normal;
out vec2 v_texcoord;
out float v_material;

uniform mat4 projMatrix;
uniform mat4 mvMatrix;
uniform mat3 normalMatrix;

void main() {
   v_vertex = o_vertex.xyz;
   v_normal = normalMatrix * o_normal;
   v_texcoord = o_texcoord;
   v_material = o_material;

   gl_Position = projMatrix * mvMatrix * o_vertex;
}
