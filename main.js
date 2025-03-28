import {Canvas, Option} from '/js/dom.js';
import {vec3, mat4} from '/js/gl-matrix-3.4.1/index.js';
import {getIDs, removeFromParent, show, hide, hidden, asyncfor} from '/js/util.js';
import {GLUtil} from '/js/gl-util.js';
import {makeGradient} from '/js/gl-util-Gradient.js';
import {Mouse3D} from '/js/mouse3d.js';
const ids = getIDs();
const urlparams = new URLSearchParams(window.location.search);

//populated in R.init
let gl;
let glutil;
let canvas;
let rotMat;
let shader;
let vtxBuf;
let mouse;

let R = new function() {
	let thiz = this;

	this.init = function() {
		
		//build dom

		canvas = Canvas({
			style : {
				zIndex : '-1',
				position : 'absolute',
			},
			appendTo : document.body,
		});

		//get gl context

		try {
			glutil = new GLUtil({canvas:canvas});
			gl = glutil.context;
		} catch (e) {
			removeFromParent(canvas);
			show(ids.webglfail);
			throw e;
		}
		glutil.import('Gradient', makeGradient);
		show(ids.menu);
		
		glutil.view.fovY = 45;
		//glutil.onfps = function(fps) { ids.fps.innerText = fps + ' fps'; };
		
		rotMat = mat4.create();
		mat4.identity(rotMat);

		glutil.view.pos[2] = 3;

		//create shaders

		shader = new glutil.Program({
			vertexCode : `
in vec3 vtx;
uniform mat4 mvMat;
uniform mat4 projMat;
out vec3 pos;
void main() {	
	pos = vtx;
	gl_Position = projMat * mvMat * vec4(vtx, 1.0); 
}
`,
			fragmentCode : `
in vec3 pos;
uniform sampler2D volTex, hsvTex;
uniform float dz;
uniform float solidThreshold;
uniform float alphaGamma;
uniform float valueAlpha;
out vec4 fragColor;
void main() { 
	float value = texture(volTex, vec2(pos.xy * .5 + .5)).r;
	fragColor = texture(hsvTex, vec2(value, .5));
	if (value > solidThreshold) {
		fragColor.w = 1.;
	} else {
		//alpha gamma: dz := dz ^ (gamma / (1+epsilon - gamma))
		fragColor.w = pow(dz, alphaGamma / (1.0001 - alphaGamma));
		//valueAlpha = 0 : alpha is constant(alphaGamma)
		//valueAlpha = 1 : alpha is a function of value such that value 0 <=> nothing, value 1 <=> original alpha
		fragColor.w *= mix(value, 1., valueAlpha);
	}
	
	//fragColor = vec4(texture(volTex, pos.xy*.5+.5).rgb, 1.);
}
`,
			uniforms : {
				dz : 1/wavefunction.dim,
				volTex : 0,
				hsvTex : 1
			}
		});
		gl.useProgram(shader.obj);

		//create buffers

		let vtxs = [
			-1,1,0,
			-1,-1,0,
			1,1,0,
			1,-1,0,
		];
		
		vtxBuf = new glutil.ArrayBuffer({
			data : vtxs
		});
		
		//init draw
		
		gl.activeTexture(gl.TEXTURE0);
		gl.clearColor(0,0,0,1);
		gl.enable(gl.BLEND);
		gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);

		//draw
		
		window.addEventListener('resize', resize);
		
		let tmpRotMat = mat4.create();	
		mouse = new Mouse3D({
			pressObj : canvas,
			move : function(dx,dy) {
				mat4.identity(tmpRotMat);
				mat4.rotate(tmpRotMat, tmpRotMat, 
					Math.PI / 180 * Math.sqrt(dx*dx + dy*dy),
					[dy, dx, 0]);
				//mat4.translate(mvMat, mvMat, [10*dx/canvas.width, -10*dy/canvas.height, 0]);
				mat4.mul(rotMat, tmpRotMat, rotMat);
				glutil.draw();
			},
			zoom : function(dz) {
				glutil.view.fovY *= Math.exp(-.0003 * dz);
				glutil.view.fovY = Math.clamp(glutil.view.fovY, 1, 179);
				glutil.updateProjection();
				glutil.draw();
			}
		});
	};
};

function factQuotient(top, bottom) {
	if (top > bottom) {
		let p = 1;
		for (let i = bottom+1; i <= top; i++) {
			p *= i;
		}
		return p;
	} else {
		let p = 1;
		for (let i = top+1; i < bottom; i++) {
			p *= i;
		}
		return 1/p;
	}
}

// assocLegendreCosThetaTable[bottomIndex][topIndex]
let associatedLegendreCosThetaTable = {
	'0' : {
		'0' : function(theta) { return 1; }
	},
	'1' : {
		'0' : function(theta) { return Math.cos(theta); },
		'1' : function(theta) { return -Math.sin(theta); }
	},
	'2' : {
		'0' : function(theta) { return .5 * (3 * Math.pow(Math.cos(theta),2) - 1); },
		'1' : function(theta) { return -3 * Math.cos(theta) * Math.sin(theta); },
		'2' : function(theta) { return 3 * Math.pow(Math.sin(theta),2); }
	},
	'3' : {
		'0' : function(theta) { return .5 * (5 * Math.pow(Math.cos(theta),2) - 3) * Math.cos(theta); },
		'1' : function(theta) { return -1.5 * (5 * Math.pow(Math.cos(theta),2) - 1) * Math.sin(theta); },
		'2' : function(theta) { return 15 * Math.cos(theta) * Math.pow(Math.sin(theta),2); },
		'3' : function(theta) { return -15 * Math.pow(Math.sin(theta),3); }
	},
	'4' : {
		'0' : function(theta) { return .125 * (35 * Math.pow(Math.cos(theta),4) - 30 * Math.pow(Math.cos(theta),2) + 3); },
		'1' : function(theta) { return -2.5 * (7 * Math.pow(Math.cos(theta),3) - 3 * Math.cos(theta)) * Math.sin(theta); },
		'2' : function(theta) { return 7.5 * (7 * Math.pow(Math.cos(theta),2) - 1) * Math.pow(Math.sin(theta),2); },
		'3' : function(theta) { return -105 * Math.cos(theta) * Math.pow(Math.sin(theta),3); },
		'4' : function(theta) { return 105 * Math.pow(Math.sin(theta),4); }
	}
};

function associatedLegendreCosTheta(l,m,theta) {
	if (m < 0) {
		m = -m
		return Math.pow(-1, m) * factQuotient(l - m, l + m) * associatedLegendreCosThetaTable[l][m](theta);
	} else {
		return associatedLegendreCosThetaTable[l][m](theta);
	}
}

//http://www.mpi-hd.mpg.de/astrophysik/HEA/internal/Numerical_Recipes/f6-8.pdf
//0 <= m <= l, -1 <= x <= 1
function associatedLegendre(l,m,x) {
	if (m < -l || m > l || Math.abs(x) > 1) throw 'bad';
	if (m < 0) {
		m = -m
		return Math.pow(-1, m) * factQuotient(l - m, l + m) * associatedLegendre(l,m,x);
	}
	
	let pmm = 1.;
	if (m > 0) {
		let somx2 = Math.sqrt((1-x)*(1-x));
		let fact = 1.;
		for (let i = 1; i <= m; i++) {
			pmm = -pmm * fact * somx2;
			fact = fact + 2;
		}
	}
	
	if (l == m) {
		return pmm;
	}
	pmmp1 = x * (2*m+1)*pmm;
	if (l == m+1) {
		return pmmp1;
	}
	for (let ll = m+2; ll >= 1; ll--) {
		pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
		pmm = pmmp1;
		pmmp1 = pll;
	}
	return pll;
}

function sphericalHarmonic(result, theta, phi, l, m) {
	let magn = Math.sqrt(
		(2*l+1) / (4*Math.PI)
		* factQuotient(l-m, l+m)
	) 
	//* associatedLegendre(l,m,Math.cos(theta));
	* associatedLegendreCosTheta(l,m,theta);
	let arg = m * phi;
	result[0] = magn * Math.cos(arg);
	result[1] = magn * Math.sin(arg);
}

function ncr(n,k) {
	let prod = 1;
	for (let i = 1; i <= k; i++) {
		prod *= (n - k + i) / i;
	}
	return prod;
}

function Laguerre(x, alpha, n) {
	let sum = 0;
	for (let i = 0; i <= n; i++) {
		sum += Math.pow(-1, i) * ncr(n + alpha, n - 1) * Math.pow(x, i) / factQuotient(i,1);
	}
	return sum;
}

function hydrogen(result, r, theta, phi, n, l, m) {
	let a0 = 1;	//.0529nm = hbar / (me * c * alpha)
	let b = 2 / (n * a0);
	let rho = b * r;
	let magn = Math.sqrt(
		b * b * b
		* factQuotient(n-l-1, n+l) / (2*n)
	) * Math.exp(-.5 * rho)
	* Math.pow(rho,l)
	* Laguerre(rho, 2*l+1, n-l-1);
	sphericalHarmonic(result, theta, phi, l, m);
	result[0] *= magn;
	result[1] *= magn;
}

let wavefunction = new function() {
	this.dim = 64;
	this.n = 4;	//orbit param
	this.l = 2;	//sh param
	this.m = 0;	//sh param
	this.min = vec3.clone([-50,-50,-50]);
	this.max = vec3.clone([50,50,50]);
	this.size = vec3.create();
	vec3.sub(this.size, this.max, this.min);
	this.init = function() {
		this.hsvTex = new glutil.Gradient.GradientTexture({
			width : 256, 
			colors : [
				[0, 0, 0],
				[0, 0, 1],
				[1, 0, 1],
				[1, 0, 0],
				[1, 1, 0],
				[1, 1, 1]
			],
			dontRepeat : true
		});

		//let data = new Uint8Array(this.dim * this.dim * 3);
		this.slices = [[],[],[]];	//for each axii 
		for (let dim = 0; dim < 3; dim++) {
			let dim1 = (dim+1)%3;
			let dim2 = (dim+2)%3;
			for (let w = 0; w < this.dim; w++) {
				let slice = {}; 
				slice.tex = new glutil.Texture2D();
				this.slices[dim].push(slice);
			}
		}
		
		this.rebuild();
	};

	this.rebuild = function() {
		show(ids.show_calculating);
		ids.calculating.setAttribute('value', '0');
		
		this.minv = Infinity;
		this.maxv = -Infinity;
		let floatData = [];//new Float32Array(this.dim * this.dim); 
		let result = [0,0];

		if (this.building !== undefined) {
			//interrupt the previous build 
			clearInterval(this.building);
		}

		let x = vec3.create();
		let thiz = this;
		this.building = asyncfor({
			start : 0,
			end : 3 * this.dim, 
			timeout : 100,
			callback : function(q) {
				let w = q % thiz.dim;
				let dim = (q - w) / thiz.dim;
				if (floatData[dim] === undefined) floatData[dim] = [];
				let dim1 = (dim+1)%3;
				let dim2 = (dim+2)%3;
			
				floatData[dim][w] = new Float32Array(thiz.dim * thiz.dim);
				let slice = thiz.slices[dim][w]; 
				for (let u = 0; u < thiz.dim; u++) {
					for (let v = 0; v < thiz.dim; v++) {
						/**/
						x[dim] = (w+.5)/thiz.dim * (thiz.max[dim] - thiz.min[dim]) + thiz.min[dim];
						x[dim1] = (u+.5)/thiz.dim * (thiz.max[dim1] - thiz.min[dim1]) + thiz.min[dim1];
						x[dim2] = (v+.5)/thiz.dim * (thiz.max[dim2] - thiz.min[dim2]) + thiz.min[dim2];
						let r = vec3.length(x);
						let phi = Math.atan2(x[1], x[0]);
						let theta = Math.acos(x[2] / r);
						//let v = sphericalHarmonic(theta, phi, thiz.l, thiz.m);
						hydrogen(result, r, theta, phi, thiz.n, thiz.l, thiz.m);
						let value = Math.sqrt(result[0]*result[0] + result[1]*result[1]);
						floatData[dim][w][u+thiz.dim*v] = value;
						if (value < thiz.minv) thiz.minv = value;
						if (value > thiz.maxv) thiz.maxv = value;
						/**/	
						/*	
						x[dim] = (w+.5)/thiz.dim * 255; 
						x[dim1] = (u+.5)/thiz.dim * 255; 
						x[dim2] = (v+.5)/thiz.dim * 255;
						data[0+3*(u+thiz.dim*v)] = x[0];
						data[1+3*(u+thiz.dim*v)] = x[1];
						data[2+3*(u+thiz.dim*v)] = x[2];
						*/
					}
				}
				ids.calculating.setAttribute('value', 100*q/(3*thiz.dim-1));
			},
			done : function() {
				let byteData = new Uint8Array(thiz.dim * thiz.dim);
				for (let dim = 0; dim < 3; dim++) {
					for (let w = 0; w < thiz.dim; w++) {
						for (let u = 0; u < thiz.dim * thiz.dim; ++u) {
							byteData[u] = parseInt((floatData[dim][w][u] - thiz.minv) / (thiz.maxv - thiz.minv) * 255);
						}
						let slice = thiz.slices[dim][w]; 
						slice.tex.bind();
						slice.tex.setArgs({
							width:thiz.dim,
							height:thiz.dim,
							format:gl.LUMINANCE,
							internalFormat:gl.LUMINANCE,
							type:gl.UNSIGNED_BYTE,
							data:byteData,
							minFilter:gl.NEAREST,
							magFilter:gl.LINEAR,
							wrap:{
								s:gl.CLAMP_TO_EDGE,
								t:gl.CLAMP_TO_EDGE
							}
						});
						slice.tex.unbind();
					}
				}
				thiz.ready = true;
				hide(ids.show_calculating);
				glutil.draw();
			}
		});
	};
};

let slideThreshold = 1.;
let alphaGamma = .25;
let valueAlpha = 0.;

function resize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	glutil.resize();
	glutil.draw();

	let info = ids.info;
	let width = window.innerWidth 
		- parseInt(info.style.paddingLeft)
		- parseInt(info.style.paddingRight);
	info.style.width = width+'px';
	let height = window.innerHeight
		- parseInt(info.style.paddingTop)
		- parseInt(info.style.paddingBottom);
	info.style.height = (height - 32)+'px';

}

function drawQuad(mv) {
	shader.setUniform('mvMat', mv);
	gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);	//4 = vtxcount
}

let viewMat = mat4.create();	
let objMat = mat4.create();

ids.panelButton.addEventListener('click', e => { 
	if (hidden(ids.panel)) {
		show(ids.panel);
		hide(ids.info);
	} else {
		hide(ids.panel);
	}
});
ids.infoButton.addEventListener('click', e => {
	if (hidden(ids.info)) {
		show(ids.info);
		hide(ids.panel);
	} else {
		hide(ids.info);
	}
});


ids.solid_threshold_slider.value = 100*(1 - slideThreshold);
ids.solid_threshold_slider.addEventListener('input', e => {
	slideThreshold = 1 - ids.solid_threshold_slider.value/100;
	glutil.draw();
});
ids.alpha_gamma_slider.value = 100*(1 - alphaGamma);
ids.alpha_gamma_slider.addEventListener('input', e => {
	alphaGamma = 1 - ids.alpha_gamma_slider.value/100;
	glutil.draw();
});
ids.value_alpha_slider.value = 100*(1 - valueAlpha);
ids.value_alpha_slider.addEventListener('input', e => {
	valueAlpha = 1 - ids.value_alpha_slider.value/100;
	glutil.draw();
});

R.init();

//N goes from 0 to 5, L goes from 0 to N-1, M goes from -L to L
const waveparams = {};
const permutations = ids.permutations;
for (let n = 0; n <= 5; ++n) {
	for (let l = 0; l < n; ++l) {
		for (let m = -l; m <= l; ++m) {
			let text = 'N='+n+' L='+l+' M='+m;
			let option = Option({innerText:text, appendTo:permutations});
			if (n == wavefunction.n && l == wavefunction.l && m == wavefunction.m) {
				option.setAttribute('selected', 'true');
			}
			waveparams[text] = {n:n, l:l, m:m};
		}
	}
}
permutations.addEventListener('change', e => {
	let params = waveparams[permutations.value];
	wavefunction.n = params.n; 
	wavefunction.l = params.l; 
	wavefunction.m = params.m; 
	wavefunction.rebuild(); 
});

glutil.ondraw = function() {
	if (!wavefunction.ready) return;
	
	shader.setAttr('vtx', vtxBuf);
	shader.setUniform('projMat', glutil.scene.projMat);

	wavefunction.hsvTex.bind(1);
	gl.activeTexture(gl.TEXTURE0);
	
	mat4.multiply(viewMat, glutil.scene.mvMat, rotMat);
	//now pick the dir and order (front vs back)
	// based on the major axis
	// (highest z-component of each axis)
	let dir = 0;
	let bestZ = viewMat[2];
	for (let i = 1; i < 3; i++) {
		let z = viewMat[4*i+2];
		if (Math.abs(z) > Math.abs(bestZ)) {
			bestZ = z;
			dir = i;
		}
	}

	let firstI, lastI, stepI;
	if (bestZ < 0) {
		firstI = wavefunction.slices[dir].length-1;
		lastI = -1;
		stepI = -1;
	} else {
		firstI = 0;
		lastI = wavefunction.slices[dir].length;
		stepI = 1;
	}

	switch (dir) {
	case 0://x-align: 
		mat4.rotate(viewMat, viewMat, Math.PI/2, [0,1,0]);
		mat4.rotate(viewMat, viewMat, Math.PI/2, [0,0,1]);
		break;	
	case 1:	//y-align:
		mat4.rotate(viewMat, viewMat, Math.PI/2, [-1,0,0]);
		mat4.rotate(viewMat, viewMat, Math.PI/2, [0,0,-1]);
		break;
	case 2:
		break;	//z-aligned is default
	}

	shader.setUniforms({
		solidThreshold : slideThreshold,
		alphaGamma : alphaGamma,
		valueAlpha : valueAlpha
	});
	for (let i = firstI; i != lastI; i+=stepI) {
		mat4.translate(objMat, viewMat, [0,0,2*(i/(wavefunction.slices[dir].length-1))-1]);
		let slice = wavefunction.slices[dir][i];
		gl.bindTexture(gl.TEXTURE_2D, slice.tex.obj);
		drawQuad(objMat);
	}
}

wavefunction.init();

resize();
