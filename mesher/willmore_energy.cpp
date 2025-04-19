#include "PeriodicMesher.h"
#include <Eigen/PardisoSupport>

using namespace msf;

using mint = int;
using mreal = double;


Eigen::SparseMatrix<double> DerivativeAssembler(PeriodSurfaceMesh& m, mreal weight)
{
	//auto fInds = mesh->getFaceIndices();
	//VertexIndices vInds = mesh->getVertexIndices();

	mint vertex_count = m.n_vertices();
	mint primitive_count = m.n_faces();
	mint primitive_length = 3;

	auto outer = std::vector<mint>(primitive_count * primitive_length + 1, 0);
	auto inner = std::vector<mint>(primitive_count * primitive_length, 0);
	auto values = std::vector<mreal>(primitive_count * primitive_length, weight);

	for (auto face : m.faces())
	{
		//mint i = fInds[face];
		mint i = face.idx();

		auto he = face.halfedge();

		mint i0 = he.from().idx();
		mint i1 = he.next().from().idx();
		mint i2 = he.next().next().from().idx();

		outer[primitive_length * i + 1] = primitive_length * i + 1;
		outer[primitive_length * i + 2] = primitive_length * i + 2;
		outer[primitive_length * i + 3] = primitive_length * i + 3;

		inner[primitive_length * i + 0] = i0;
		inner[primitive_length * i + 1] = i1;
		inner[primitive_length * i + 2] = i2;
	}
	return Eigen::Map<Eigen::SparseMatrix<mreal>>(vertex_count, primitive_count * primitive_length, primitive_count * primitive_length, &outer[0], &inner[0], &values[0]);
}


void WillmoreEnergy_Differential(PeriodSurfaceMesh& m, Eigen::MatrixXd& output_diff, double& output_value)
{
	auto [x, primitives] = m.getVFlist();
	//auto primitives = getPrimitiveIndices(mesh, geom);
	//auto x = getVertexPositions(mesh, geom);

	mint vertex_count = x.rows();
	mint dim = x.cols();
	mint primitive_count = primitives.rows();
	mint primitive_length = primitives.cols();
	//requireMeanCurvatureVectors();

	auto buffer = Eigen::MatrixXd(primitive_count * primitive_length, dim);

	Eigen::MatrixX3d H(m.n_vertices(), 3); H.setZero();
	Eigen::VectorXd H_squared(m.n_vertices()); H_squared.setZero();

	output_diff.resize(m.n_vertices(), 3); output_diff.setZero();
	output_value = 0;

	for (auto vh : m.vertices()) {
		auto [o, ring] = m.find1ring(vh, 1);
		auto [Hv, A] = m.meanH(o, ring);
		//Hv /= -2 * A;
		Hv /= 2 * A;
		H.row(vh.idx()) = Hv.transpose();
		H_squared[vh.idx()] = Hv.squaredNorm();
		output_value += H_squared[vh.idx()] * A;
	}

	mreal weight = 1;

	mreal one_third = 1. / 3.;

#pragma omp parallel
	for (mint i = 0; i < primitive_count; ++i)
	{
		mint i0 = primitives(i, 0);
		mint i1 = primitives(i, 1);
		mint i2 = primitives(i, 2);

		Eigen::Vector3d x0 = x.row(i0).transpose();
		Eigen::Vector3d x1 = x.row(i1).transpose();
		Eigen::Vector3d x2 = x.row(i2).transpose();
		x1 = x0 + make_period((x1 - x0).eval(), 2, 1);
		x2 = x0 + make_period((x2 - x0).eval(), 2, 1);

		mreal x00 = x0[0];
		mreal x01 = x0[1];
		mreal x02 = x0[2];

		mreal x10 = x1[0];
		mreal x11 = x1[1];
		mreal x12 = x1[2];

		mreal x20 = x2[0];
		mreal x21 = x2[1];
		mreal x22 = x2[2];

		mreal v00 = H(i0, 0);
		mreal v01 = H(i0, 1);
		mreal v02 = H(i0, 2);

		mreal v10 = H(i1, 0);
		mreal v11 = H(i1, 1);
		mreal v12 = H(i1, 2);

		mreal v20 = H(i2, 0);
		mreal v21 = H(i2, 1);
		mreal v22 = H(i2, 2);

		mreal weight = one_third * (H_squared(i0) + H_squared(i1) + H_squared(i2));

		mreal s0 = -(x01 * x10);
		mreal s1 = x00 * x11;
		mreal s2 = x01 * x20;
		mreal s3 = -(x11 * x20);
		mreal s4 = -(x00 * x21);
		mreal s5 = x10 * x21;
		mreal s6 = s0 + s1 + s2 + s3 + s4 + s5;
		mreal s7 = s6 * s6;
		mreal s8 = x02 * x10;
		mreal s9 = -(x00 * x12);
		mreal s10 = -(x02 * x20);
		mreal s11 = x12 * x20;
		mreal s12 = x00 * x22;
		mreal s13 = -(x10 * x22);
		mreal s14 = s10 + s11 + s12 + s13 + s8 + s9;
		mreal s15 = s14 * s14;
		mreal s16 = -(x02 * x11);
		mreal s17 = x01 * x12;
		mreal s18 = x02 * x21;
		mreal s19 = -(x12 * x21);
		mreal s20 = -(x01 * x22);
		mreal s21 = x11 * x22;
		mreal s22 = s16 + s17 + s18 + s19 + s20 + s21;
		mreal s23 = s22 * s22;
		mreal s24 = s15 + s23 + s7;
		mreal s25 = sqrt(s24);
		mreal s26 = 1 / s25;
		mreal s27 = -x21;
		mreal s28 = s27 + x11;
		mreal s29 = 2 * s28 * s6;
		mreal s30 = -x12;
		mreal s31 = s30 + x22;
		mreal s32 = 2 * s14 * s31;
		mreal s33 = s29 + s32;
		mreal s34 = s24 * s25;
		mreal s35 = 1 / s34;
		mreal s36 = -x10;
		mreal s37 = s36 + x20;
		mreal s38 = -x00;
		mreal s39 = s38 + x10;
		mreal s40 = -x20;
		mreal s41 = s40 + x00;
		mreal s42 = s40 + x10;
		mreal s43 = -x11;
		mreal s44 = -x02;
		mreal s45 = s43 + x01;
		mreal s46 = s44 + x12;
		mreal s47 = -x22;
		mreal s48 = -x01;
		mreal s49 = s48 + x21;
		mreal s50 = s47 + x02;
		mreal s51 = s33 * s33;
		mreal s52 = s28 * s28;
		mreal s53 = s31 * s31;
		mreal s54 = s38 + x20;
		mreal s55 = s36 + x00;
		mreal s56 = 2 * s37 * s6;
		mreal s57 = s47 + x12;
		mreal s58 = 2 * s22 * s57;
		mreal s59 = s56 + s58;
		mreal s60 = -(s33 * s35 * s59) / 8.;
		mreal s61 = (s26 * s28 * s37) / 2.;
		mreal s62 = s60 + s61;
		mreal s63 = 2 * s49 * s6;
		mreal s64 = 2 * s14 * s50;
		mreal s65 = s63 + s64;
		mreal s66 = -2 * s6;
		mreal s67 = 2 * s45 * s6;
		mreal s68 = 2 * s14 * s46;
		mreal s69 = s67 + s68;
		mreal s70 = 2 * s6;
		mreal s71 = 2 * s39 * s6;
		mreal s72 = s30 + x02;
		mreal s73 = 2 * s22 * s72;
		mreal s74 = s71 + s73;
		mreal s75 = s59 * s59;
		mreal s76 = s37 * s37;
		mreal s77 = s57 * s57;
		mreal s78 = 2 * s14 * s42;
		mreal s79 = s43 + x21;
		mreal s80 = 2 * s22 * s79;
		mreal s81 = s78 + s80;
		mreal s82 = 2 * s41 * s6;
		mreal s83 = s44 + x22;
		mreal s84 = 2 * s22 * s83;
		mreal s85 = s82 + s84;
		mreal s86 = 2 * s14 * s55;
		mreal s87 = s48 + x11;
		mreal s88 = 2 * s22 * s87;
		mreal s89 = s86 + s88;
		mreal s90 = 2 * s14 * s54;
		mreal s91 = s27 + x01;
		mreal s92 = 2 * s22 * s91;
		mreal s93 = s90 + s92;
		mreal s94 = s81 * s81;
		mreal s95 = s42 * s42;
		mreal s96 = s79 * s79;
		mreal s97 = -(s35 * s59 * s81) / 8.;
		mreal s98 = (s26 * s57 * s79) / 2.;
		mreal s99 = s97 + s98;
		mreal s100 = -(s33 * s35 * s81) / 8.;
		mreal s101 = (s26 * s31 * s42) / 2.;
		mreal s102 = s100 + s101;
		mreal s103 = -2 * s14;
		mreal s104 = 2 * s14;
		mreal s105 = -2 * s22;
		mreal s106 = 2 * s22;
		mreal s107 = -(s35 * s59 * s65) / 8.;
		mreal s108 = 2 * s37 * s49;
		mreal s109 = s108 + s66;
		mreal s110 = (s109 * s26) / 4.;
		mreal s111 = s107 + s110;
		mreal s112 = s65 * s65;
		mreal s113 = s49 * s49;
		mreal s114 = s50 * s50;
		mreal s115 = -(s33 * s35 * s65) / 8.;
		mreal s116 = 2 * s28 * s49;
		mreal s117 = 2 * s31 * s50;
		mreal s118 = s116 + s117;
		mreal s119 = (s118 * s26) / 4.;
		mreal s120 = s115 + s119;
		mreal s121 = -(s35 * s65 * s81) / 8.;
		mreal s122 = 2 * s42 * s50;
		mreal s123 = s104 + s122;
		mreal s124 = (s123 * s26) / 4.;
		mreal s125 = s121 + s124;
		mreal s126 = -(s35 * s65 * s85) / 8.;
		mreal s127 = (s26 * s41 * s49) / 2.;
		mreal s128 = s126 + s127;
		mreal s129 = -(s33 * s35 * s85) / 8.;
		mreal s130 = 2 * s28 * s41;
		mreal s131 = s130 + s70;
		mreal s132 = (s131 * s26) / 4.;
		mreal s133 = s129 + s132;
		mreal s134 = -(s35 * s59 * s85) / 8.;
		mreal s135 = 2 * s37 * s41;
		mreal s136 = 2 * s57 * s83;
		mreal s137 = s135 + s136;
		mreal s138 = (s137 * s26) / 4.;
		mreal s139 = s134 + s138;
		mreal s140 = s85 * s85;
		mreal s141 = s41 * s41;
		mreal s142 = s83 * s83;
		mreal s143 = -(s35 * s81 * s85) / 8.;
		mreal s144 = 2 * s79 * s83;
		mreal s145 = s105 + s144;
		mreal s146 = (s145 * s26) / 4.;
		mreal s147 = s143 + s146;
		mreal s148 = s93 * s93;
		mreal s149 = s54 * s54;
		mreal s150 = s91 * s91;
		mreal s151 = -(s35 * s81 * s93) / 8.;
		mreal s152 = 2 * s42 * s54;
		mreal s153 = 2 * s79 * s91;
		mreal s154 = s152 + s153;
		mreal s155 = (s154 * s26) / 4.;
		mreal s156 = s151 + s155;
		mreal s157 = -(s35 * s65 * s93) / 8.;
		mreal s158 = (s26 * s50 * s54) / 2.;
		mreal s159 = s157 + s158;
		mreal s160 = -(s35 * s85 * s93) / 8.;
		mreal s161 = (s26 * s83 * s91) / 2.;
		mreal s162 = s160 + s161;
		mreal s163 = -(s33 * s35 * s93) / 8.;
		mreal s164 = 2 * s31 * s54;
		mreal s165 = s103 + s164;
		mreal s166 = (s165 * s26) / 4.;
		mreal s167 = s163 + s166;
		mreal s168 = -(s35 * s59 * s93) / 8.;
		mreal s169 = 2 * s57 * s91;
		mreal s170 = s106 + s169;
		mreal s171 = (s170 * s26) / 4.;
		mreal s172 = s168 + s171;
		mreal s173 = s69 * s69;
		mreal s174 = s45 * s45;
		mreal s175 = s46 * s46;
		mreal s176 = -(s35 * s69 * s85) / 8.;
		mreal s177 = 2 * s41 * s45;
		mreal s178 = s177 + s66;
		mreal s179 = (s178 * s26) / 4.;
		mreal s180 = s176 + s179;
		mreal s181 = -(s35 * s59 * s69) / 8.;
		mreal s182 = 2 * s37 * s45;
		mreal s183 = s182 + s70;
		mreal s184 = (s183 * s26) / 4.;
		mreal s185 = s181 + s184;
		mreal s186 = -(s35 * s65 * s69) / 8.;
		mreal s187 = 2 * s45 * s49;
		mreal s188 = 2 * s46 * s50;
		mreal s189 = s187 + s188;
		mreal s190 = (s189 * s26) / 4.;
		mreal s191 = s186 + s190;
		mreal s192 = -(s33 * s35 * s69) / 8.;
		mreal s193 = 2 * s28 * s45;
		mreal s194 = 2 * s31 * s46;
		mreal s195 = s193 + s194;
		mreal s196 = (s195 * s26) / 4.;
		mreal s197 = s192 + s196;
		mreal s198 = -(s35 * s69 * s81) / 8.;
		mreal s199 = 2 * s42 * s46;
		mreal s200 = s103 + s199;
		mreal s201 = (s200 * s26) / 4.;
		mreal s202 = s198 + s201;
		mreal s203 = -(s35 * s69 * s93) / 8.;
		mreal s204 = 2 * s46 * s54;
		mreal s205 = s104 + s204;
		mreal s206 = (s205 * s26) / 4.;
		mreal s207 = s203 + s206;
		mreal s208 = -(s35 * s69 * s74) / 8.;
		mreal s209 = (s26 * s39 * s45) / 2.;
		mreal s210 = s208 + s209;
		mreal s211 = s74 * s74;
		mreal s212 = s39 * s39;
		mreal s213 = s72 * s72;
		mreal s214 = -(s33 * s35 * s74) / 8.;
		mreal s215 = 2 * s28 * s39;
		mreal s216 = s215 + s66;
		mreal s217 = (s216 * s26) / 4.;
		mreal s218 = s214 + s217;
		mreal s219 = -(s35 * s65 * s74) / 8.;
		mreal s220 = 2 * s39 * s49;
		mreal s221 = s220 + s70;
		mreal s222 = (s221 * s26) / 4.;
		mreal s223 = s219 + s222;
		mreal s224 = -(s35 * s59 * s74) / 8.;
		mreal s225 = 2 * s37 * s39;
		mreal s226 = 2 * s57 * s72;
		mreal s227 = s225 + s226;
		mreal s228 = (s227 * s26) / 4.;
		mreal s229 = s224 + s228;
		mreal s230 = -(s35 * s74 * s85) / 8.;
		mreal s231 = 2 * s39 * s41;
		mreal s232 = 2 * s72 * s83;
		mreal s233 = s231 + s232;
		mreal s234 = (s233 * s26) / 4.;
		mreal s235 = s230 + s234;
		mreal s236 = -(s35 * s74 * s93) / 8.;
		mreal s237 = 2 * s72 * s91;
		mreal s238 = s105 + s237;
		mreal s239 = (s238 * s26) / 4.;
		mreal s240 = s236 + s239;
		mreal s241 = -(s35 * s74 * s81) / 8.;
		mreal s242 = 2 * s72 * s79;
		mreal s243 = s106 + s242;
		mreal s244 = (s243 * s26) / 4.;
		mreal s245 = s241 + s244;
		mreal s246 = s89 * s89;
		mreal s247 = s55 * s55;
		mreal s248 = s87 * s87;
		mreal s249 = -(s35 * s74 * s89) / 8.;
		mreal s250 = (s26 * s72 * s87) / 2.;
		mreal s251 = s249 + s250;
		mreal s252 = -(s35 * s69 * s89) / 8.;
		mreal s253 = (s26 * s46 * s55) / 2.;
		mreal s254 = s252 + s253;
		mreal s255 = -(s35 * s89 * s93) / 8.;
		mreal s256 = 2 * s54 * s55;
		mreal s257 = 2 * s87 * s91;
		mreal s258 = s256 + s257;
		mreal s259 = (s258 * s26) / 4.;
		mreal s260 = s255 + s259;
		mreal s261 = -(s35 * s81 * s89) / 8.;
		mreal s262 = 2 * s42 * s55;
		mreal s263 = 2 * s79 * s87;
		mreal s264 = s262 + s263;
		mreal s265 = (s26 * s264) / 4.;
		mreal s266 = s261 + s265;
		mreal s267 = -(s35 * s65 * s89) / 8.;
		mreal s268 = 2 * s50 * s55;
		mreal s269 = s103 + s268;
		mreal s270 = (s26 * s269) / 4.;
		mreal s271 = s267 + s270;
		mreal s272 = -(s33 * s35 * s89) / 8.;
		mreal s273 = 2 * s31 * s55;
		mreal s274 = s104 + s273;
		mreal s275 = (s26 * s274) / 4.;
		mreal s276 = s272 + s275;
		mreal s277 = -(s35 * s59 * s89) / 8.;
		mreal s278 = 2 * s57 * s87;
		mreal s279 = s105 + s278;
		mreal s280 = (s26 * s279) / 4.;
		mreal s281 = s277 + s280;
		mreal s282 = -(s35 * s85 * s89) / 8.;
		mreal s283 = 2 * s83 * s87;
		mreal s284 = s106 + s283;
		mreal s285 = (s26 * s284) / 4.;
		mreal s286 = s282 + s285;
		buffer(3 * i + 0, 0) = -((-(s35 * s51) / 8. + (s26 * (2 * s52 + 2 * s53)) / 4.) * v00) - s62 * v01 - s102 * v02 - s120 * v10 - s133 * v11 - s167 * v12 - s197 * v20 - s218 * v21 - s276 * v22 - (s26 * s33 * weight) / 4.;
		buffer(3 * i + 0, 1) = -(s62 * v00) - (-(s35 * s75) / 8. + (s26 * (2 * s76 + 2 * s77)) / 4.) * v01 - s99 * v02 - s111 * v10 - s139 * v11 - s172 * v12 - s185 * v20 - s229 * v21 - s281 * v22 - (s26 * s59 * weight) / 4.;
		buffer(3 * i + 0, 2) = -(s102 * v00) - s99 * v01 - (-(s35 * s94) / 8. + (s26 * (2 * s95 + 2 * s96)) / 4.) * v02 - s125 * v10 - s147 * v11 - s156 * v12 - s202 * v20 - s245 * v21 - s266 * v22 - (s26 * s81 * weight) / 4.;
		buffer(3 * i + 1, 0) = -(s120 * v00) - s111 * v01 - s125 * v02 - (((2 * s113 + 2 * s114) * s26) / 4. - (s112 * s35) / 8.) * v10 - s128 * v11 - s159 * v12 - s191 * v20 - s223 * v21 - s271 * v22 - (s26 * s65 * weight) / 4.;
		buffer(3 * i + 1, 1) = -(s133 * v00) - s139 * v01 - s147 * v02 - s128 * v10 - (((2 * s141 + 2 * s142) * s26) / 4. - (s140 * s35) / 8.) * v11 - s162 * v12 - s180 * v20 - s235 * v21 - s286 * v22 - (s26 * s85 * weight) / 4.;
		buffer(3 * i + 1, 2) = -(s167 * v00) - s172 * v01 - s156 * v02 - s159 * v10 - s162 * v11 - (((2 * s149 + 2 * s150) * s26) / 4. - (s148 * s35) / 8.) * v12 - s207 * v20 - s240 * v21 - s260 * v22 - (s26 * s93 * weight) / 4.;
		buffer(3 * i + 2, 0) = -(s197 * v00) - s185 * v01 - s202 * v02 - s191 * v10 - s180 * v11 - s207 * v12 - (((2 * s174 + 2 * s175) * s26) / 4. - (s173 * s35) / 8.) * v20 - s210 * v21 - s254 * v22 - (s26 * s69 * weight) / 4.;
		buffer(3 * i + 2, 1) = -(s218 * v00) - s229 * v01 - s245 * v02 - s223 * v10 - s235 * v11 - s240 * v12 - s210 * v20 - (((2 * s212 + 2 * s213) * s26) / 4. - (s211 * s35) / 8.) * v21 - s251 * v22 - (s26 * s74 * weight) / 4.;
		buffer(3 * i + 2, 2) = -(s276 * v00) - s281 * v01 - s266 * v02 - s271 * v10 - s286 * v11 - s260 * v12 - s254 * v20 - s251 * v21 - (((2 * s247 + 2 * s248) * s26) / 4. - (s246 * s35) / 8.) * v22 - (s26 * s89 * weight) / 4.;
	}

	output_diff += weight * DerivativeAssembler(m, 1) * buffer;
}

std::tuple<double, Eigen::MatrixXd> msf::PeriodSurfaceMesh::willmore_energy(void)
{
	Eigen::MatrixXd diff;
	double value;
	WillmoreEnergy_Differential(*this, diff, value);
	return { value, diff };
}

extern std::tuple<Eigen::VectorXd, std::vector<Eigen::Vector3d>> eval_vertex_mass(PeriodSurfaceMesh& m);

extern std::vector<double> aux_number;

void test_willmore_energy(void) {
#if 0
	std::string mfile = "D:/projects/minisurf/image/siggraph/adc-opt/P0.5-31.stl";
	PeriodSurfaceMesh m;
	m.readMergePeriodBoundary(mfile);
#else
	std::string mfile = "D:/projects/minisurf/image/siggraph/adc-opt/ellip123.obj";
	PeriodSurfaceMesh m;
	m.read(mfile, false, false, false);
#endif
	m._nvdof = m.n_vertices();
	double eng_last = 1e30;
	double t_step = aux_number[3];
	for (int iter = 0; iter < aux_number[4]; iter++) {
		std::cout << "*iter " << iter << std::endl;
		if (iter % 50 == 0) {
			//m.save("iter.obj");
			m.savePeriodicMesh("iter.obj", std::set<OM::SmartVertexHandle>{}, 1);
		}
		auto [value, diff] = m.willmore_energy();
		m._nvdof = m.n_vertices();
		if (aux_number[0] == 1) {
			auto [Av, Nv] = eval_vertex_mass(m);
			Eigen::SparseMatrix<double> L = -0.2 * m.getPeriodicLaplacian(2, 1);
			L += Av.asDiagonal();
			Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(L);
			diff = so.solve((Av.asDiagonal() * diff).eval());
			Eigen::Vector3d u_mean = (Av.asDiagonal() * diff).colwise().sum().transpose() / Av.sum();
			diff.rowwise() -= u_mean.transpose();
		}
		else if (aux_number[0] == 2) {
			auto [Av, Nv] = eval_vertex_mass(m);
			Eigen::SparseMatrix<double> L = m.getPeriodicLaplacian(2, 1);
			Eigen::SparseMatrix<double> L2 = L * Av.cwiseInverse().asDiagonal() * L;
			Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so;
			if (aux_number[5] == 1) {
				so.compute(L);  
			} else if (aux_number[5] == 2) {
				so.compute(L2); // invalid result!
			}
			//Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(L2);
			diff = so.solve((Av.asDiagonal() * diff).eval());
			Eigen::Vector3d u_mean = (Av.asDiagonal() * diff).colwise().sum().transpose() / Av.sum();
			diff.rowwise() -= u_mean.transpose();
		}
		std::cout << "diff size = " << diff.rows() << ", " << diff.cols() << std::endl;
		std::cout << "willmore energy = " << value << std::endl;
		if (value > eng_last || iter % 50 == 0) {
			m.save("befrm.obj");
			std::cout << "Remeshing...";
			m.delaunayRemesh(10, 0.03, 0.02);
			std::cout << "Finished" << std::endl;
			m.save("aftrm.obj");
			eng_last = 1e30;
			continue;
		}
		eng_last = value;
		{
			//std::ofstream ofs("willdiff", std::ios::binary);
			//for (auto vh : m.vertices()) {
			//	auto p = m.point(vh);
			//	ofs.write((const char*)p.data(), sizeof(p));
			//	Eigen::Vector3d v = diff.row(vh.idx()).transpose();
			//	ofs.write((const char*)v.data(), sizeof(v));
			//}
		}
		for (auto vh : m.vertices()) {
			//std::cout << diff.row(vh.idx()) << std::endl;
			//std::cout << toOM((diff.row(vh.idx()) * t_step).eval()) << std::endl;
			m.point(vh) -= toOM((diff.row(vh.idx()) * t_step).eval());
		}
	}
}

Eigen::MatrixX3d willmore_H2_gradient(msf::PeriodSurfaceMesh& m, const Eigen::SparseMatrix<double>& L, Eigen::VectorXd Av, std::vector<Eigen::Vector3d>& Nv)
{
	auto [value, diff] = m.willmore_energy();
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(L);
	diff = so.solve((Av.asDiagonal() * diff).eval());
	Eigen::Vector3d u_mean = (Av.asDiagonal() * diff).colwise().sum().transpose() / Av.sum();
	diff.rowwise() -= u_mean.transpose();
	return diff;
}


Eigen::MatrixX3d willmore_H2_gradient(msf::PeriodSurfaceMesh& m)
{
	auto [Av, Nv] = eval_vertex_mass(m);
	Eigen::SparseMatrix<double> L = m.getPeriodicLaplacian(2, 1);
	return willmore_H2_gradient(m, L, Av, Nv);
}
