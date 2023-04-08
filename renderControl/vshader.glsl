// origin data
attribute vec4 o_vertex;
attribute vec3 o_normal;
attribute vec2 o_texcoord;
attribute float o_material;

// vshader data
varying vec3 v_vertex;
varying vec3 v_normal;
varying vec2 v_texcoord;
varying float v_material;

// uniform data
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
