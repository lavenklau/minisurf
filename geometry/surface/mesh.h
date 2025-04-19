#pragma once

#include "surface.h"

#define _USE_MATH_DEFINES
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include<OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
//#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include<vector>
#include "facet/BezierTri.h"


namespace OpenMesh {
	struct MeshTraits : public OpenMesh::DefaultTraits
	{
		typedef OpenMesh::Vec3d Point;
		typedef OpenMesh::Vec3d Normal;
		typedef OpenMesh::Vec2d TexCoord2D;

		VertexAttributes(OpenMesh::Attributes::Status);
		FaceAttributes(OpenMesh::Attributes::Status);
		EdgeAttributes(OpenMesh::Attributes::Status);
		HalfedgeAttributes(OpenMesh::Attributes::Status);
	};

}


BEGIN_MINSURF_NAMESPACE

namespace OM = ::OpenMesh;

typedef OpenMesh::TriMesh_ArrayKernelT<OpenMesh::MeshTraits> OM_Mesh;

template<typename Raw>
struct SmartH { using type = Raw; };
template<> struct SmartH<OM::VertexHandle> { using type = OM::SmartVertexHandle; };
template<> struct SmartH<OM::HalfedgeHandle> { using type = OM::SmartHalfedgeHandle; };
template<> struct SmartH<OM::FaceHandle> { using type = OM::SmartFaceHandle; };


template<typename P> Eigen::Vector3<Real> toEigen(const P& p) { return Eigen::Vector3<Real>(p[0], p[1], p[2]); }
template<typename P> auto toOM(const P& p) { return OpenMesh::Vec3d(p[0], p[1], p[2]); }

class MeshSurface 
	: public OM_Mesh
{
public:
	virtual std::pair<Eigen::Vector3<Real>, Eigen::Vector3<Real>> boundingBox(void) const {
		Eigen::Vector3<Real> pmin, pmax;
		pmin.setConstant(1e30);
		pmax.setConstant(-1e30);
		for (auto vh : vertices()) {
			auto vp = point(vh);
			pmax = pmax.cwiseMax(Eigen::Vector3<Real>(vp[0], vp[1], vp[2]));
			pmin = pmin.cwiseMin(Eigen::Vector3<Real>(vp[0], vp[1], vp[2]));
		}
		return { pmin, pmax };
	}
	Eigen::Matrix<Real, 3, 3> getFaceVertex(OM::FaceHandle fh) const {
		Eigen::Matrix3<Real> m;
		int counter = 0;
		for (auto fv : fv_ccw_range(fh)) {
			auto p = point(fv);
			for (int i = 0; i < 3; i++) {
				m(i, counter) = p[i];
			}
			counter++;
		}
		return m;
	}
	std::array<OM::SmartVertexHandle, 3> getFaceVertexHandle(OM::FaceHandle fh) const {
		std::array<OM::SmartVertexHandle, 3> vhlist;
		int counter = 0;
		for (auto fv : fv_ccw_range(fh)) {
			if (!fv.is_valid()) {
				printf("\033[31m Invalid face vertex\033[0m\n");
			}
			vhlist[counter++] = fv;
		}
		//if (counter < 3) {
		//	printf("\033[31m Warning : incomplete face\033[0m\n");
		//}
		return vhlist;
	}
	std::array<OM::SmartHalfedgeHandle, 3> getFaceHeHandle(OM::FaceHandle fh) const {
		std::array<OM::SmartHalfedgeHandle, 3> helist;
		int counter = 0;
		for (auto fhe : fh_range(fh)) {
			helist[counter++] = fhe;
		}
		return helist;
	}

	virtual std::pair<Eigen::Matrix<msf::Real, -1, 3>, Eigen::Matrix<int, -1, 3>> getVFlist(void);

	msf::OM::SmartHalfedgeHandle findHe(OM::VertexHandle from, OM::VertexHandle to);

	int heid(OM::HalfedgeHandle h) const;

	void request_normal_status(void);

	Real cothe(OM::HalfedgeHandle he) const;

	Real cote(OM::EdgeHandle he) const;

	virtual bool read(std::string filepath, bool normalize);
	virtual void normalize(void);
	virtual bool read(const Eigen::MatrixX3<Real>& V, const Eigen::MatrixX3i& F);
	void write(std::string filename);
	Real totalArea(void) const;
	std::vector<OM::FaceHandle> facelist(void) const;
	std::vector<OM::VertexHandle> vertexlist(void) const;

	void faceVhandle(OM::FaceHandle fh, OM::VertexHandle vh[3]) const;

	std::pair<OM::SmartHalfedgeHandle, OM::SmartHalfedgeHandle> adjacent_he(OM::FaceHandle f1, OM::FaceHandle f2) const;
	
	std::vector<std::vector<OM::SmartHalfedgeHandle>> extract_boundary_loop(void) const;

	Real average_edge_length(void);
};

class BezierMesh
	: public MeshSurface 
{
protected:
	// size = [Face, 3] 
	Eigen::MatrixX3i matching;
	Eigen::MatrixX3<Real> ctrl_points;
	Eigen::MatrixXi f_cid;
	int order;
	int _zid(Real u, Real v) const;
	int _zid(int i, int j, int k) const;
	//std::vector<VertexFlag> ctrl_pflags;
protected:
	// edge order clock wise, [v=0] -> [u=0] -> [1-u-v=0] 
	Eigen::Vector2<Real> transition(OM::HalfedgeHandle hefrom, OM::HalfedgeHandle heto, Real u, Real v) const;
	Eigen::Vector3i transition(OM::HalfedgeHandle hefrom, OM::HalfedgeHandle heto, int i, int j, int k) const;
	int start_vid(OM::HalfedgeHandle h, bool global = false);
	int end_vid(OM::HalfedgeHandle h, bool global = false);
public:
	int heid(OM::HalfedgeHandle h) const;
	auto& getControlPoints(void) const { return ctrl_points; }
	auto& getFaceCPlist(void) const { return f_cid; }
	void updateControlPoints(const Eigen::MatrixX3<Real>& pmat) { ctrl_points = pmat; }
	void set_order(int order_);
	int n_ctrlpoints(void);
	// compute control points 
	virtual void buildBezierMesh(void);

	BezierTri<Real, -1> getBezierTri(OM::FaceHandle fh);
};

END_MINSURF_NAMESPACE
