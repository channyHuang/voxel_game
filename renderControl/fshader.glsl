varying highp vec3 v_vertex;
varying highp vec3 v_normal;
varying vec2 v_texcoord;
varying float v_material;

uniform highp vec3 lightPos;
//uniform sampler2D texture;
uniform sampler2D textures[17];
//uniform sampler3D textures;

void main() {
   highp vec3 L = normalize(lightPos - v_vertex);
   highp float NL = max(dot(normalize(v_normal), L), 0.0);

//   highp vec3 color = vec3(0.39, 1.0, 0.0);
//   highp vec3 col = clamp(color * 0.8 + color * 0.2 * NL, 0.0, 1.0);
//   gl_FragColor = vec4(col, 1.0);

//   highp vec4 tex = texture2D(textures[0], v_texcoord);
//   gl_FragColor = clamp(tex * 0.8 + tex * 0.2 * NL, 0.0, 1.0);

   int i;
   int material = int(v_material);
   int matCount = 0;
   highp vec4 mixTexture = vec4(0.0, 0.0, 0.0, 0.0);
   for (i = 0; i < 17; ++i) {
       if (material & (1 << i)) {
            highp vec4 tex = texture2D(textures[i], v_texcoord);
            mixTexture += tex;
            matCount++;
       }
   }
   mixTexture /= float(matCount);

   gl_FragColor = clamp(mixTexture * 0.8 + mixTexture * 0.2 * NL, 0.0, 1.0);
}
