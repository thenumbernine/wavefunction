//populated in R.init
var canvas;
var rotMat;
var shader;
var vtxBuf;
var gl;
var mouse;

var R = new function() {
	var thiz = this;

	this.init = function() {
		
		//build dom

		canvas = $('<canvas>').appendTo(document.body).get(0);

		//get gl context

		try {
			gl = GL.init(canvas);
		} catch (e) {
			$(canvas).remove();
			$('#webglfail').show();
			throw e;
		}
		$('#menu').show();

		GL.view.fovY = 45;
		//GL.onfps = function(fps) { $('#fps').text(fps + ' fps'); };
		
		rotMat = mat4.create();
		mat4.identity(rotMat);

		GL.view.pos[2] = 3;

		//create shaders

		shader = new GL.ShaderProgram({
			vertexPrecision : 'best',
			vertexCode : mlstr(function(){/*
attribute vec3 vtx;
uniform mat4 mvMat;
uniform mat4 projMat;
varying vec3 pos;
void main() {	
	pos = vtx;
	gl_Position = projMat * mvMat * vec4(vtx, 1.0); 
}
*/}),
			fragmentPrecision : 'best',
			fragmentCode : mlstr(function(){/*
varying vec3 pos;
uniform sampler2D volTex, hsvTex;
uniform float dz;
uniform float solidThreshold;
uniform float alphaGamma;
uniform float valueAlpha;
void main() { 
	float value = texture2D(volTex, vec2(pos.xy * .5 + .5)).r;
	gl_FragColor = texture2D(hsvTex, vec2(value, .5));
	if (value > solidThreshold) {
		gl_FragColor.w = 1.;
	} else {
		//alpha gamma: dz := dz ^ (gamma / (1+epsilon - gamma))
		gl_FragColor.w = pow(dz, alphaGamma / (1.0001 - alphaGamma));
		//valueAlpha = 0 : alpha is constant(alphaGamma)
		//valueAlpha = 1 : alpha is a function of value such that value 0 <=> nothing, value 1 <=> original alpha
		gl_FragColor.w *= mix(value, 1., valueAlpha);
	}
	
	//gl_FragColor = vec4(texture2D(volTex, pos.xy*.5+.5).rgb, 1.);
}
*/}),
			uniforms : {
				dz : 1/wavefunction.dim,
				volTex : 0,
				hsvTex : 1
			}
		});
		gl.useProgram(shader.obj);

		//create buffers

		var vtxs = [
			-1,1,0,
			-1,-1,0,
			1,1,0,
			1,-1,0,
		];
		
		vtxBuf = new GL.ArrayBuffer({data:vtxs});
		
		//init draw
		
		gl.activeTexture(gl.TEXTURE0);
		gl.clearColor(0,0,0,1);
		gl.enable(gl.BLEND);
		gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);

		//draw
		
		$(window).resize(resize);
		
		var tmpRotMat = mat4.create();	
		mouse = new Mouse3D({
			pressObj : canvas,
			move : function(dx,dy) {
				mat4.identity(tmpRotMat);
				mat4.rotate(tmpRotMat, tmpRotMat, 
					Math.PI / 180 * Math.sqrt(dx*dx + dy*dy),
					[dy, dx, 0]);
				//mat4.translate(mvMat, mvMat, [10*dx/canvas.width, -10*dy/canvas.height, 0]);
				mat4.mul(rotMat, tmpRotMat, rotMat);
				GL.draw();
			},
			zoom : function(dz) {
				GL.view.fovY *= Math.exp(-.0003 * dz);
				GL.view.fovY = Math.clamp(GL.view.fovY, 1, 179);
				GL.updateProjection();
				GL.draw();
			}
		});
	};
};

function factQuotient(top, bottom) {
	if (top > bottom) {
		var p = 1;
		for (var i = bottom+1; i <= top; i++) {
			p *= i;
		}
		return p;
	} else {
		var p = 1;
		for (var i = top+1; i < bottom; i++) {
			p *= i;
		}
		return 1/p;
	}
}

// assocLegendreCosThetaTable[bottomIndex][topIndex]
var associatedLegendreCosThetaTable = {
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
	
	var pmm = 1.;
	if (m > 0) {
		var somx2 = Math.sqrt((1-x)*(1-x));
		var fact = 1.;
		for (var i = 1; i <= m; i++) {
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
	for (var ll = m+2; ll >= 1; ll--) {
		pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
		pmm = pmmp1;
		pmmp1 = pll;
	}
	return pll;
}

function sphericalHarmonic(result, theta, phi, l, m) {
	var magn = Math.sqrt(
		(2*l+1) / (4*Math.PI)
		* factQuotient(l-m, l+m)
	) 
	//* associatedLegendre(l,m,Math.cos(theta));
	* associatedLegendreCosTheta(l,m,theta);
	var arg = m * phi;
	result[0] = magn * Math.cos(arg);
	result[1] = magn * Math.sin(arg);
}

function ncr(n,k) {
	var prod = 1;
	for (var i = 1; i <= k; i++) {
		prod *= (n - k + i) / i;
	}
	return prod;
}

function Laguerre(x, alpha, n) {
	var sum = 0;
	for (var i = 0; i <= n; i++) {
		sum += Math.pow(-1, i) * ncr(n + alpha, n - 1) * Math.pow(x, i) / factQuotient(i,1);
	}
	return sum;
}

function hydrogen(result, r, theta, phi, n, l, m) {
	var a0 = 1;	//.0529nm = hbar / (me * c * alpha)
	var b = 2 / (n * a0);
	var rho = b * r;
	var magn = Math.sqrt(
		b * b * b
		* factQuotient(n-l-1, n+l) / (2*n)
	) * Math.exp(-.5 * rho)
	* Math.pow(rho,l)
	* Laguerre(rho, 2*l+1, n-l-1);
	sphericalHarmonic(result, theta, phi, l, m);
	result[0] *= magn;
	result[1] *= magn;
}

var wavefunction = new function() {
	this.dim = 256;
	this.n = 4;	//orbit param
	this.l = 2;	//sh param
	this.m = 0;	//sh param
	this.min = vec3.clone([-50,-50,-50]);
	this.max = vec3.clone([50,50,50]);
	this.size = vec3.create();
	vec3.sub(this.size, this.max, this.min);
	this.init = function() {
		this.hsvTex = new GL.GradientTexture({
			width:256, 
			colors:[
				[0,0,0],
				[0,0,1],
				[1,0,1],
				[1,0,0],
				[1,1,0],
				[1,1,1]
			],
			dontRepeat : true
		});

		//var data = new Uint8Array(this.dim * this.dim * 3);
		this.slices = [[],[],[]];	//for each axii 
		for (var dim = 0; dim < 3; dim++) {
			var dim1 = (dim+1)%3;
			var dim2 = (dim+2)%3;
			for (var w = 0; w < this.dim; w++) {
				var slice = {}; 
				slice.tex = new GL.Texture2D();
				this.slices[dim].push(slice);
			}
		}
		
		this.rebuild();
	};

	this.rebuild = function() {
		$('#show-calculating').show();
		$('#calculating').attr('value', '0');
		
		this.minv = Infinity;
		this.maxv = -Infinity;
		var floatData = [];//new Float32Array(this.dim * this.dim); 
		var result = [0,0];

		if (this.building !== undefined) {
			//interrupt the previous build 
			clearInterval(this.building);
		}

		var x = vec3.create();
		var thiz = this;
		this.building = asyncfor({
			start : 0,
			end : 3 * this.dim, 
			timeout : 100,
			callback : function(q) {
				var w = q % thiz.dim;
				var dim = (q - w) / thiz.dim;
				if (floatData[dim] === undefined) floatData[dim] = [];
				var dim1 = (dim+1)%3;
				var dim2 = (dim+2)%3;
			
				floatData[dim][w] = new Float32Array(thiz.dim * thiz.dim);
				var slice = thiz.slices[dim][w]; 
				for (var u = 0; u < thiz.dim; u++) {
					for (var v = 0; v < thiz.dim; v++) {
						/**/
						x[dim] = (w+.5)/thiz.dim * (thiz.max[dim] - thiz.min[dim]) + thiz.min[dim];
						x[dim1] = (u+.5)/thiz.dim * (thiz.max[dim1] - thiz.min[dim1]) + thiz.min[dim1];
						x[dim2] = (v+.5)/thiz.dim * (thiz.max[dim2] - thiz.min[dim2]) + thiz.min[dim2];
						var r = vec3.length(x);
						var phi = Math.atan2(x[1], x[0]);
						var theta = Math.acos(x[2] / r);
						//var v = sphericalHarmonic(theta, phi, thiz.l, thiz.m);
						hydrogen(result, r, theta, phi, thiz.n, thiz.l, thiz.m);
						var value = Math.sqrt(result[0]*result[0] + result[1]*result[1]);
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
				$('#calculating').attr('value', 100*q/(3*thiz.dim-1));
			},
			done : function() {
				var byteData = new Uint8Array(thiz.dim * thiz.dim);
				for (var dim = 0; dim < 3; dim++) {
					for (var w = 0; w < thiz.dim; w++) {
						for (var u = 0; u < thiz.dim * thiz.dim; ++u) {
							byteData[u] = parseInt((floatData[dim][w][u] - thiz.minv) / (thiz.maxv - thiz.minv) * 255);
						}
						var slice = thiz.slices[dim][w]; 
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
				$('#show-calculating').hide();
				GL.draw();
			}
		});
	};
};

var slideThreshold = 1.;
var alphaGamma = .25;
var valueAlpha = 0.;
$(document).ready(function() {
	$('#panelButton').click(function() {
		var panel = $('#panel');	
		if (panel.css('display') == 'none') {
			panel.show();
			$('#info').hide();
		} else {
			panel.hide();
		}
	});
	$('#infoButton').click(function() {
		var info = $('#info');
		if (info.css('display') == 'none') {
			info.show();
			$('#panel').hide();
		} else {
			info.hide();
		}
	});
	
	
	$('#solid-threshold-slider').slider({
		range : 'max',
		width : '200px',
		min : 0,
		max : 100,
		value : 100*(1 - slideThreshold),
		slide : function(event, ui) {
			slideThreshold = 1 - ui.value/100;
			GL.draw();
		}
	});
	$('#alpha-gamma-slider').slider({
		range : 'max',
		width : '200px',
		min : 0,
		max : 100,
		value : 100*(1 - alphaGamma),
		slide : function(event, ui) {
			alphaGamma = 1 - ui.value/100;
			GL.draw();
		}
	});
	$('#value-alpha-slider').slider({
		range : 'max',
		width : '200px',
		min : 0,
		max : 100,
		value : 100*(1 - valueAlpha),
		slide : function(event, ui) {
			valueAlpha = 1 - ui.value/100;
			GL.draw();
		}
	});

	R.init();

	//N goes from 0 to 5, L goes from 0 to N-1, M goes from -L to L
	var waveparams = {};
	var permutations = $('#permutations');
	for (var n = 0; n <= 5; ++n) {
		for (var l = 0; l < n; ++l) {
			for (var m = -l; m <= l; ++m) {
				var text = 'N='+n+' L='+l+' M='+m;
				var option = $('<option>', {text:text});
				option.appendTo(permutations);
				if (n == wavefunction.n && l == wavefunction.l && m == wavefunction.m) {
					option.attr('selected', 'true');
				}
				waveparams[text] = {n:n, l:l, m:m};
			}
		}
	}
	permutations.change(function() {
		var params = waveparams[permutations.val()];
		wavefunction.n = params.n; 
		wavefunction.l = params.l; 
		wavefunction.m = params.m; 
		wavefunction.rebuild(); 
	});

	GL.ondraw = function() {
		if (!wavefunction.ready) return;
		
		shader.setAttr('vtx', vtxBuf);
		shader.setUniform('projMat', GL.projMat);

		wavefunction.hsvTex.bind(1);
		gl.activeTexture(gl.TEXTURE0);
		
		mat4.multiply(viewMat, GL.mvMat, rotMat);
		//now pick the dir and order (front vs back)
		// based on the major axis
		// (highest z-component of each axis)
		var dir = 0;
		var bestZ = viewMat[2];
		for (var i = 1; i < 3; i++) {
			var z = viewMat[4*i+2];
			if (Math.abs(z) > Math.abs(bestZ)) {
				bestZ = z;
				dir = i;
			}
		}

		var firstI, lastI, stepI;
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
		for (var i = firstI; i != lastI; i+=stepI) {
			mat4.translate(objMat, viewMat, [0,0,2*(i/(wavefunction.slices[dir].length-1))-1]);
			var slice = wavefunction.slices[dir][i];
			gl.bindTexture(gl.TEXTURE_2D, slice.tex.obj);
			drawQuad(objMat);
		}
	}
	
	wavefunction.init();

	resize();
})

function resize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	GL.resize();
	GL.draw();

	var info = $('#info');
	var width = window.innerWidth 
		- parseInt(info.css('padding-left'))
		- parseInt(info.css('padding-right'));
	info.width(width);
	var height = window.innerHeight
		- parseInt(info.css('padding-top'))
		- parseInt(info.css('padding-bottom'));
	info.height(height - 32);

}

function drawQuad(mv) {
	shader.setUniform('mvMat', mv);
	gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);	//4 = vtxcount
}

var viewMat = mat4.create();	
var objMat = mat4.create();

