varying highp vec3 vert;
varying highp vec3 vertNormal;
varying vec2 v_texcoord;

uniform highp vec3 lightPos;
uniform sampler2D texture;
uniform sampler2D texture1;
//uniform sampler3D textures;

void main() {
   highp vec3 L = normalize(lightPos - vert);
   highp float NL = max(dot(normalize(vertNormal), L), 0.0);

   //highp vec3 color = vec3(0.39, 1.0, 0.0);
   //highp vec3 col = clamp(color * 0.2 + color * 0.8 * NL, 0.0, 1.0);
   //gl_FragColor = vec4(col, 1.0);

   highp vec4 tex = texture2D(texture, v_texcoord);
   highp vec4 tex1 = mix(texture2D(texture, v_texcoord), texture2D(texture1, v_texcoord), 0.5);
   //gl_FragColor = clamp(tex * 0.2 + tex * 0.8 * NL, 0.0, 1.0);
   gl_FragColor = clamp(tex1 * 0.8 + tex1 * 0.2 * NL, 0.0, 1.0);
}
