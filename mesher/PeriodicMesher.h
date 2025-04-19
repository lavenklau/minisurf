#pragma once
#include "Config.h"
#include "geometry/surface/mesh.h"
#include <unordered_map>
#include "bbox/BBox.h"
#include "Flag.h"
#include "grid/PeriodicGrid.h"

BEGIN_MINSURF_NAMESPACE

class PeriodSurfaceMesh
	: public MeshSurface 
{

public:
	struct HETag {
		bool pass_boundary[3];
		Eigen::Vector3<Real> torus_trans;
		HETag(void) { for (int i = 0; i < 3; i++) { pass_boundary[i] = false; torus_trans[i] = 0; } }
		void updateFlag(void) { for (int i = 0; i < 3; i++) { pass_boundary[i] = torus_trans[i] != 0; } }
		HETag(const Eigen::Vector3<Real>& trans) : torus_trans(trans) { for (int i = 0; i < 3; i++) { pass_boundary[i] = trans[i] != 0; } }
		bool has_trans(void) const { return pass_boundary[0] + pass_boundary[1] + pass_boundary[2]; }
		HETag opp(void) const { return HETag(-torus_trans); }
	};

	// ** ** ** **** **** **** **** **** ** TO BE DEPRECATED * * * * * * * ** * * ** * * ** * * ** * * ** * * ** * * *
	// the second is reduced dof
	std::unordered_map<OpenMesh::Vec3d, int> _vpid;
	std::unordered_map<OpenMesh::VertexHandle, int> _vdofid;
	// period vertices
	std::vector<std::vector<OpenMesh::VertexHandle>> _vPeriodGroup;
	// halfedge transfer across periodic boundary
	std::unordered_map<OpenMesh::SmartHalfedgeHandle, OpenMesh::SmartHalfedgeHandle> _heTransfer;
	// unreduced vector with period dof
	std::vector<VertexFlag> _vflags;
	// periodic edge to adjacent faces
	std::unordered_map<OM::EdgeHandle, std::pair<OM::FaceHandle, OM::FaceHandle>> _egperface;
	int _nvdof;
	// bounding box of original model
	BBox _bbox;
	// for periodic search
	//std::unique_ptr<aabb::Tree> _tree;
	// ** ** ** **** **** **** **** **** ** * * * * * * * * ** * * ** * * ** * * ** * * ** * * ** * * *
private:
	void setPeriod(void);
	void normalizeBb(void);
	std::vector<Eigen::Vector3<Real>> scaleTo(const BBox& new_bb);
	void scaleTo(const BBox& new_bb, std::vector<Eigen::Vector3<Real>>& new_vlist);
	void scaleTo(const BBox& new_bb, Eigen::MatrixX3<Real>& new_vlist);
	void updateBb(void);
	//OM::SmartHalfedgeHandle findHe(OM::VertexHandle from, OM::VertexHandle to);
	void relaxInTorusH1(const std::vector<std::vector<size_t>>& cutlocus, const std::vector<int>& handle_tag);
	void relaxInTorusH2(const std::vector<std::vector<size_t>>& cutlocus, const std::vector<int>& handle_tag);
	void relaxInTorusARAP(const std::vector<std::vector<size_t>>& cutlocus, const std::vector<int>& handle_tag);
	void meanCurvatureFlow(Eigen::MatrixX3<Real>& vlist, const std::vector<OM::SmartVertexHandle>& vcut, const std::vector<OM::SmartHalfedgeHandle>& hetrans, const std::vector<HETag>& hetags);
	Eigen::Matrix3<msf::Real> matchRot(const std::vector<Real>& wij, const Eigen::Vector3<Real>& o, const std::vector<Eigen::Vector3<Real>>& ring, const Eigen::Vector3<Real>& o1, const std::vector<Eigen::Vector3<Real>>& ring1);
	Eigen::Matrix3<msf::Real> matchRot(const Eigen::Vector3<Real>& o, const std::vector<Eigen::Vector3<Real>>& ring, const Eigen::Vector3<Real>& o1, const std::vector<Eigen::Vector3<Real>>& ring1);
	void relaxInTorusLaplacian(const std::vector<std::vector<size_t>>& cutlocus, const std::vector<int>& handle_tag);
	void relaxInTorus(const std::vector<std::vector<size_t>>& cutlocus, const std::vector<int>& handle_tag, int type);
	void initAABBTree(void);
	Real _cot_gradient(const Eigen::Vector3<Real>& a, const Eigen::Vector3<Real>& b, Eigen::Vector3<Real>& agrad, Eigen::Vector3<Real>& bgrad) const;
	std::tuple<Real, Real> meanH2(Eigen::Vector3<Real>& o, std::vector<Eigen::Vector3<Real>>& ring, bool area_averaged);
	void meanH(const Eigen::Vector3<Real>& o, const std::vector<Eigen::Vector3<Real>>& ring, Eigen::Vector3<Real>& Ho);
	std::tuple<Real, Real> meanHs(Eigen::Vector3<Real>& o, std::vector<Eigen::Vector3<Real>>& ring, Real s, bool area_averaged);
	std::tuple<Real, Real> meanH2_test(Eigen::Vector3<Real>& o, std::vector<Eigen::Vector3<Real>>& ring, bool area_averaged);
	std::tuple<Real, Real> meanH2_nograd(Eigen::Vector3<Real>& o, std::vector<Eigen::Vector3<Real>>& ring, bool area_averaged);
	std::set<OM::SmartHalfedgeHandle> getBoundaryHalfedges(void);
	std::tuple<std::vector<OM::SmartHalfedgeHandle>, std::vector<HETag>> toTransHETag(const std::vector<std::vector<size_t>>& cutlocus, const std::vector<int>& handle_tag, Eigen::MatrixX3<Real>& vtrans);
public:
	std::tuple<Eigen::Vector3d, Real> meanH(const Eigen::Vector3<Real>& o, const std::vector<Eigen::Vector3<Real>>& ring);
	Eigen::Vector3<Real> eval_period_edge(OM::SmartHalfedgeHandle he, double t) const;
	std::tuple<Eigen::Vector3<Real>, std::vector<Eigen::Vector3<Real>>> find1ring(const Eigen::VectorX<Real>& vlist, OM::SmartVertexHandle vh, const std::vector<HETag>& hetags);
	std::tuple<Eigen::Vector3<Real>, std::vector<Eigen::Vector3<Real>>> find1ring(OM::SmartVertexHandle vh, const std::vector<HETag>& hetags);
	std::tuple<Eigen::Vector3<Real>, std::vector<Eigen::Vector3<Real>>> find1ring(OM::SmartVertexHandle vh, Real deduce) const;
	std::tuple<Eigen::Vector3<Real>, std::vector<Eigen::Vector3<Real>>> find1ring(const Eigen::MatrixX3<Real>& vlist, OM::SmartVertexHandle vh, const std::vector<HETag>& hetags);
	std::tuple<Eigen::Vector3<Real>, std::vector<Eigen::Vector3<Real>>> find1ring(const Eigen::MatrixX3<Real>& vlist, OM::SmartVertexHandle vh) const;
	std::tuple<std::vector<OM::SmartVertexHandle>,std::vector<OM::SmartHalfedgeHandle>> find1ring(OM::SmartVertexHandle vh);
	std::vector<Real> find1ring(OM::SmartVertexHandle vh, const Eigen::SparseMatrix<Real>& L);
	Eigen::Vector3<msf::Real> bary1ring(const Eigen::Vector3<Real>& v, const std::vector<Eigen::Vector3<Real>>& ring);
	// type -1 : centroid  ; -2 : mixed 
	std::tuple<Real, Eigen::Vector3<Real>> area1ring(const Eigen::Vector3<Real>& v, const std::vector<Eigen::Vector3<Real>>& ring, int type = 1);
	// d0->d1, left, right
	std::tuple<Eigen::Vector3<Real>, Eigen::Vector3<Real>, Eigen::Vector3<Real>, Eigen::Vector3<Real>>
		find1diag(OM::SmartHalfedgeHandle he, const std::vector<HETag>& hetags);
	std::tuple<Eigen::Vector3<Real>, Eigen::Vector3<Real>, Eigen::Vector3<Real>, Eigen::Vector3<Real>>
		find1diag(OM::SmartHalfedgeHandle he, Real deduce);
	std::vector<Real> cot1ring(const Eigen::Vector3<Real>& vlist, const std::vector<Eigen::Vector3<Real>>& ring);
	//std::tuple<Eigen::Vector3<Real>, Eigen::Vector3<Real>> curvature1ring(const Eigen::Vector3<Real>& v0, const std::vector<Eigen::Vector3<Real>>& ring);
	//std::tuple<Real, Eigen::Vector3<Real>, Real, Eigen::Vector3<Real>> curvaturekring(OM::SmartVertexHandle vh, int k, Real deduce = 0.5);
	std::tuple<Real, Real> normalDistribution(const Eigen::Vector3<Real>& q);
	//std::tuple<Eigen::MatrixX3<Real>, std::map<int, int>> genSmoothField(bool period, Real smoothness, Real alignment);
	Real angle(const Eigen::Vector3<Real>& va, const Eigen::Vector3<Real>& vo, const Eigen::Vector3<Real>& vb);
	Real cosangle(const Eigen::Vector3<Real>& vo, OM::SmartHalfedgeHandle he);
	//void getVVDirac(Eigen::SparseMatrix<Real>& Dv, const Eigen::MatrixX3<Real>& vlist, const Eigen::VectorX<Real> Hv, const std::vector<OM::SmartVertexHandle>& vcut, const std::vector<OM::SmartHalfedgeHandle>& hetrans, const std::vector<HETag>& hetags);
	//Eigen::VectorX<Real> MobiusBalance(const Eigen::MatrixX3<Real>& vlist, const Eigen::VectorX<Real>& vA, const Eigen::VectorX<Real>& vA_goal, const Eigen::MatrixX3<Real>& vN);
	std::vector<int> getVertexValence(void);
	std::vector<Real> getFaceArea(void);
	Eigen::Vector3<Real> getFaceArea(OM::FaceHandle fh, Real deduce);
	Eigen::Matrix3<Real> getFacePeriodVertex(OM::FaceHandle fh, Real deduce) const;
	Eigen::Matrix3<Real> getFacePeriodVertex(OM::FaceHandle fh, OM::VertexHandle vhsync, Real deduce) const;
	Eigen::VectorX<Real> getFaceArea(const Eigen::MatrixX3<Real>& vlist, const std::vector<HETag>& het);
	Real periodicArea(Real deduce);
	// type = 1 average by area, 2 average by angle
	Eigen::MatrixX3<Real> getVertexNormals(int type, const Eigen::MatrixX3<Real>& vlist, const std::vector<OM::SmartVertexHandle>& vcut, const std::vector<OM::SmartHalfedgeHandle>& hetrans, const std::vector<HETag>& hetags);
	std::tuple<Eigen::Vector3<Real>, Eigen::Vector3<Real>> getVertexNormal(const Eigen::Vector3<Real>& o, const std::vector<Eigen::Vector3<Real>>& ring, int type = 1);
	std::tuple<Eigen::MatrixX3<Real>, Eigen::MatrixX3<Real>> getFaceNormal(Real period, Real deduce);
	Eigen::MatrixX3<Real> getVertexNormal(Real period, Real deduce);
	Eigen::MatrixX3<Real> getFaceFrame(Real period, Real deduce);
	//Eigen::MatrixX3<Real> spinXform(const Eigen::SparseMatrix<Real>& L, const Eigen::MatrixX3<Real>& vold, const std::vector<Eigen::Quaternion<Real>>& v_quat,
		//const std::vector<OM::SmartVertexHandle>& vcut, const std::vector<OM::SmartHalfedgeHandle>& hetrans, const std::vector<HETag>& hetags);
	//void getFVDirac(Eigen::SparseMatrix<Real>& Df, Eigen::MatrixX3<Real>& vlist, const std::vector<OM::SmartVertexHandle>& vcut, const std::vector<OM::SmartHalfedgeHandle>& hetrans, const std::vector<HETag>& hetags);
	template<typename Bmat, typename Sol>
	Eigen::ComputationInfo solveLU(const Eigen::SparseMatrix<Real>& A, const Bmat& bmat, Sol& sol) {
		Eigen::SparseLU<Eigen::SparseMatrix<Real>> solver(A);
		sol = solver.solve(bmat);
		return solver.info();
	}
public:
	void pertubMesh(int max_iter, Real strength, unsigned int seed);
	bool hasSliver(Real angleThres = 30, Real period = 2, Real deduce = 0.5);
	Real aspect_ratio(OM::FaceHandle fh, Real period = 2, Real deduce = 0.5);
	Eigen::Matrix3<Real> face_edge_vector(OM::FaceHandle fh, Real period = 2, Real deduce = 0.5);
	template<typename Raw> auto toSmart(const Raw& H) { return typename SmartH<Raw>::type(H.idx(), this); }
	void mergePeriodEdges(void);
	std::vector<OM::SmartVertexHandle> mergePeriodVertices(void);
	void remove_dup_vertices(Real eps = 1e-5);
	static std::tuple<Eigen::MatrixX3<Real>, Eigen::MatrixX3i> remove_dup_vertices(const Eigen::MatrixX3<Real>& V, const Eigen::MatrixX3i& F, Real eps = 1e-5);
	Eigen::MatrixX3<Real> laplacianTranslationForce(const Eigen::SparseMatrix<Real>& L, const std::vector<OM::SmartHalfedgeHandle>& hetrans, const std::vector<HETag>& hetags);
	Real period_cote(OM::EdgeHandle e);
	//void willmore_flow(Eigen::MatrixX3<Real>& vlist, const std::vector<OM::SmartVertexHandle>& vcut, const std::vector<OM::SmartHalfedgeHandle>& hetrans, const std::vector<HETag>& hetags);
	//void test_willmore_flow(std::string filename);
	void cotweight(const Eigen::MatrixX3<Real>& V, Eigen::SparseMatrix<Real>& L, const std::vector<OM::SmartHalfedgeHandle>& hetrans, const std::vector<HETag>& hetags);
	//Eigen::VectorX<Real> vertexArea(Eigen::MatrixX3<Real>& vlist, const std::vector<OM::SmartVertexHandle>& vcut, const std::vector<OM::SmartHalfedgeHandle>& hetrans, const std::vector<HETag>& hetags);
	std::tuple<std::array<OM::SmartVertexHandle, 3>, std::array<Eigen::Vector3<Real>, 3>> findFaceVInc(const Eigen::MatrixX3<Real>& vlist, OM::SmartFaceHandle fh, const std::vector<HETag>& hetags);
	void savePeriodicMesh(std::string filename, const Eigen::MatrixX<Real>& V, const std::vector<OM::SmartVertexHandle>& vcut, const std::vector<HETag>& hetags);
	void savePeriodicMesh(std::string filename, const std::vector<OM::SmartVertexHandle>& vcut, Real detach = 0.5);
	void saveUnitCell(std::string filename);
	// [newv, newf, new2old]
	std::tuple<Eigen::MatrixX3<Real>, Eigen::MatrixX3i, Eigen::VectorXi> dupPeriodFaces(const Eigen::MatrixX3<Real>& vlist, const Eigen::MatrixX3i& flist, Real detach = 0.5);
	static PeriodicGridIndex generateIndexer(const Eigen::MatrixX3<Real> vlisst);
	// edge rotation
	//std::tuple<Eigen::VectorXi> combCrossField()
	Eigen::VectorXi findUniqueCorrespond(const Eigen::MatrixX3<Real>& p_uniq, const Eigen::MatrixX3<Real>& pdup);
	void savePeriodicMesh(std::string filename, const std::set<OM::SmartVertexHandle>& vcut, Real detach = 0.5);
	std::tuple<Eigen::MatrixX3<Real>, Eigen::MatrixX3i> getPeriodicMesh(void);
	static void savePeriodicMesh(std::string filename, const Eigen::MatrixX3<Real>& vlist, const Eigen::MatrixX3i& flist, Real period, Real detach /*= 0.5*/);
	void savePeriodicMesh(std::string filename, const Eigen::MatrixX<Real>& V, const std::set<OM::SmartVertexHandle>& vcut, const std::vector<HETag>& hetags);
	void sampleMeshes(int number, int reso, int sample_type, Real randomness, int remeshIter, int smooth_iter, Real edge_lengh);
	OM::HalfedgeHandle closedOpposite(OM::SmartHalfedgeHandle he);
	double edgeCot(OM::SmartHalfedgeHandle he, int vdof[4], Eigen::Vector3<Real> dofGradient[4]);
	Real edgeCot(OM::SmartHalfedgeHandle he, Real period, Real deduce);
	OM::FaceHandle oppositePerFaces(OM::EdgeHandle eh, OM::FaceHandle fh) const;
	Eigen::VectorX<Real> getVertexDofPosition(void) const;
	void setVertexDofPosition(const Eigen::VectorX<Real>& newV);
	void projectDirection(Eigen::VectorX<Real>& d);
	Real searchStepMeanCurvature(const Eigen::VectorX<Real>& u, const Eigen::VectorX<Real>& d);
	Eigen::VectorX<Real> meanCurvature(const Eigen::VectorX<Real>& v);
	Eigen::SparseMatrix<Real> getLaplacian(void);
	Eigen::SparseMatrix<Real> getPeriodicLaplacian(Real period, Real deduce);
	std::vector<Real> edgeCotList(void);
	bool read(std::string filepath, bool normalize = true, bool set_period = true, bool removeDup = true);
	bool read(const Eigen::MatrixX3<Real>& V, const Eigen::MatrixX3i& F, bool normalize = true, bool set_period = true, bool removeDup = true);
	std::vector<OM::SmartVertexHandle> readMergePeriodBoundary(std::string fielanme);
	std::vector<OM::SmartVertexHandle> readMergePeriodBoundary(std::string fielanme, bool normalize);
	std::vector<OM::SmartVertexHandle> mergePeriodBoundary(void);
	std::vector<OM::SmartVertexHandle> periodic_remesh(int max_iter, const std::vector<OM::SmartVertexHandle>& vcut, Real tgt_len, int smoot_iter, int pert_iter = 0, const MeshSurface* tgtmesh = nullptr);
	std::vector<OM::SmartVertexHandle> periodic_remesh(int max_iter, std::string filename, Real tgt_len, int smoot_iter, int pter_iter = 0);
	static std::tuple<Eigen::MatrixX3<Real>, Eigen::MatrixX3i, Eigen::VectorXi> weldMesh(const Eigen::MatrixX3<Real>& V1, const Eigen::MatrixX3i& F1, const Eigen::MatrixX3<Real>& V2, const Eigen::MatrixX3i& F2);
	bool is_cut_boundary(OM::SmartHalfedgeHandle he, const std::set<OM::SmartVertexHandle>& vcut, bool strictOnSamePlane = false);
	void periodic_smooth(int max_iter, int type, bool clamp_period_boundary = true);
	bool has_boundary(void);
	std::tuple<Eigen::MatrixX3<Real>, Eigen::VectorX<Real>>  periodic_mean_curvature(void);
	void deduce_hetag(std::set<OM::SmartVertexHandle>& vcut, std::vector<HETag>& hetags, Real detach_distance = 0.5);
	VertexFlag& getVFlag(OM::SmartVertexHandle vh);
	bool save(std::string filepath);
	int n_vdof(void) const { return _nvdof; }
	int vDofid(OM::VertexHandle vh) const { return _vdofid.at(vh); }
	void faceVdof(OM::FaceHandle fh, int gid[3]) const;
	void saveVdofVector(std::string filename, const Eigen::VectorXd& u, Eigen::Vector3i beg = Eigen::Vector3i(0, 1, 2), Eigen::Vector3i stride = Eigen::Vector3i(6, 6, 6));
	void makeMinimal(void);

	//std::pair<Eigen::Matrix<Real, -1, 3>, Eigen::Matrix<int, -1, 3>> getVFlist(void);
	Eigen::Matrix<int, -1, 3> getFlist(void);

	void embedMesh(std::string filename, std::vector<int> basepoints);

	void extendPeriod(double ratio);

	void extractAsymptoticLines(const std::vector<Eigen::Vector3d>& seeds, std::vector<std::vector<Eigen::Vector3d>>& plylines, std::vector<int>& faceid);

	void extractHarmonicField(std::vector<Eigen::VectorX<Real>>& val_vertex, std::vector<Eigen::VectorX<Real>>& vec_face);

	//Eigen::MatrixX<Real> extractHarmonicVectorBasis(const Eigen::SparseMatrix<Real>& L);

	//Eigen::MatrixX<std::complex<Real>> extractHolomorphicForm(Real period, Real deduce);

	//Eigen::MatrixX2<Real> holomorphicTexture(const Eigen::VectorX<std::complex<Real>>& form, const std::vector<std::vector<int>>& cutlocus);

	//Eigen::VectorX<std::complex<Real>> abelJacobMap(const Eigen::MatrixX<std::complex<Real>>& gforms, const std::vector<std::pair<int, int>>& divisors, const std::vector<std::vector<int>>& cutlocus);

	//std::vector<std::vector<int>> getCut2disk(void);

	//Eigen::MatrixX<std::complex<Real>> getJacobian(const Eigen::MatrixX<std::complex<Real>>& holoform);

	//Eigen::MatrixX3<Real> formOnFace(const Eigen::VectorX<Real>& w);

	//// [vid, order]
	//std::vector<std::pair<int, int>> findHoloformZero(const Eigen::VectorX<std::complex<Real>>& form);

	//Eigen::MatrixX3<std::complex<Real>> faceForm(const Eigen::VectorX<std::complex<Real>>& form);

	//Eigen::MatrixX<Real> getEdge2VertexTrans(const Eigen::MatrixX3<Real>& vlist, const Eigen::MatrixX<Real>& hvec, const std::vector<OM::SmartVertexHandle>& vcut, const std::vector<OM::SmartHalfedgeHandle>& hetrans, const std::vector<HETag>& hetags);

	//Eigen::MatrixX<Real> getEdge2FaceTrans(const Eigen::MatrixX3<Real>& vlist, const Eigen::MatrixX<Real>& hvec, const std::vector<OM::SmartVertexHandle>& vcut, const std::vector<OM::SmartHalfedgeHandle>& hetrans, const std::vector<HETag>& hetags);

	Eigen::MatrixX<Real> vecs2quats(const Eigen::MatrixX<Real>& v3list);

	Eigen::MatrixX<Real> quats2vecs(const Eigen::MatrixX<Real>& q4list);

	Eigen::VectorX<Real> assembleVertexFromDof(const Eigen::VectorX<Real>& vdof);

	std::vector<Eigen::VectorXd> inextensionalDisplacement(void);

	// 3F x 3
	Eigen::MatrixX3d second_fundamental_form(int type);

	void delaunayRemesh(int inner_iter, double tgtlen = -1, double min_len = -1, double adaptive = 1e-6, int outer_iter = 1);

	bool surgery();

	double periodEdgeLength(OM::EdgeHandle eh) const;

	double period_sector_angle(OM::HalfedgeHandle he) const;

	double period_face_area(OM::FaceHandle fh) const;

	double period_vertex_dual_area(OM::VertexHandle vh) const;

	void period_shift(int pm);

	void period_shift(void);

	void split_unit_cell(std::vector<Eigen::Vector3d>& vlist);

	void split_unit_cell(std::vector<double>& vlist);

	void split_unit_cell(void);

	std::tuple<double, Eigen::MatrixXd> willmore_energy(void);

	void clamp_period_boundary(double band = 2e-5);
};

class PeriodicMesher
{
public:
	
};

template<typename T>
struct Vector3Hash {
	std::size_t operator()(const Eigen::Vector3<T>& v) const {
		std::size_t hx = std::hash<T>()(v.x());
		std::size_t hy = std::hash<T>()(v.y());
		std::size_t hz = std::hash<T>()(v.z());
		return hx ^ (hy << 1) ^ (hz << 2); // Combine hashes
	}
};

template<typename T>
struct Vector3PairHash {
	std::size_t operator()(const std::pair<Eigen::Vector3<T>, Eigen::Vector3<T>>& p) const {
		std::size_t h1 = Vector3Hash<T>()(p.first);
		std::size_t h2 = Vector3Hash<T>()(p.second);
		return h1 ^ (h2 << 1); // Combine the two hashes
	}
};

// Macro for checking NaN in an Eigen matrix or vector
#define CHECK_NAN(matrix) \
    do { \
        if ((matrix).hasNaN()) { \
            std::cerr << "Warning: NaN found in " #matrix " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        } \
    } while (0)

END_MINSURF_NAMESPACE
