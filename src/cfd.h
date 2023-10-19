#ifndef THESIS_V3_CFD_H
#define THESIS_V3_CFD_H

#include "definitions.h"

CGAL::Bbox_3 computational_domain(double wCell);

CGAL::Bbox_3 Refinement_box_1(double wCell);
CGAL::Bbox_3 Refinement_box_2(double wCell);
CGAL::Bbox_3 Refinement_box_3(double wCell);
void RefinementBoxes(double wCell);

vector<int> Buildings_in_RoI(bool target, double x, double y);
Point3 centroid_building(int id); //??

void BuildingVolume(bool target, double x, double y);
void BuildingSeparation(bool target, double x, double y, double threshold);
void PedestrianCells(bool target, double x, double y, double threshold, int it);

void Meshing(bool target, double x, double y, double threshold);
void Meshing3(bool target, double x, double y, double threshold);

void write_blockMeshDict();
void write_snappyHexMeshDict(string filename);
void write_snappyHexMeshDict3(string filename);
#endif
