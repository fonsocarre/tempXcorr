//
//  Frame.cpp
//  Stitch
//
//  Created by Alfonso del Carre on 29/05/2015.
//  Copyright (c) 2015 Alfonso del Carre. All rights reserved.
//

#include "Frame.h"

PIV::Frame::Frame()
{
    //std::cerr << "Don't use this constructor!" << std::endl;
    //std::cerr << "From PIV::Frame() @ line " << __LINE__ << std::endl;
}

PIV::Frame::Frame(const H5File& file, std::string dataset_name,
                  int iFrame)
{
    this->iFrame = iFrame;

    Group picture = file.openGroup(dataset_name);
    
    Attribute att = picture.openAttribute("nx");
    DataType ntype = PredType::NATIVE_INT;
    att.read(ntype, &(this->nx));
    
    att = picture.openAttribute("ny");
    ntype = PredType::NATIVE_INT;
    att.read(ntype, &(this->ny));
    
    DataSet dataset;
    DataSpace dataspace;
    
    this->dataset_name = dataset_name;
    //x
    dataset = file.openDataSet(dataset_name + "x");
    dataspace = dataset.getSpace();
    //this->type = this->dataset.getTypeClass();
    auto type = PredType::NATIVE_DOUBLE;
    dataspace.getSimpleExtentDims(this->dims, NULL);
    this->x.resize(this->dims[0]);
    // TODO write type inference instead of just the macro
    dataset.read(&(this->x[0]), type, dataspace);
    
    //y
    dataset = file.openDataSet(dataset_name + "y");
    dataspace = dataset.getSpace();
    //this->type = this->dataset.getTypeClass();
    dataspace.getSimpleExtentDims(this->dims, NULL);
    this->y.resize(this->dims[0]);
    // TODO write type inference instead of just the macro
    dataset.read(&(this->y[0]), type, dataspace);
    
    //vx
    dataset = file.openDataSet(dataset_name + "vx");
    dataspace = dataset.getSpace();
    //this->type = this->dataset.getTypeClass();
    dataspace.getSimpleExtentDims(this->dims, NULL);
    this->vx.resize(this->dims[0]);
    // TODO write type inference instead of just the macro
    dataset.read(&(this->vx[0]), type, dataspace);
    
    //vy
    dataset = file.openDataSet(dataset_name + "vy");
    dataspace = dataset.getSpace();
    //this->type = this->dataset.getTypeClass();
    dataspace.getSimpleExtentDims(this->dims, NULL);
    this->vy.resize(this->dims[0]);
    // TODO write type inference instead of just the macro
    dataset.read(&(this->vy[0]), type, dataspace);
    
    this->n = static_cast<int>(this->dims[0]);
}

PIV::Frame& PIV::Frame::operator=(PIV::Frame &&other)
{
    this->x = std::move(other.x);
    this->y = std::move(other.y);
    this->vx = std::move(other.vx);
    this->vy = std::move(other.vy);
    
    this->dataset_name = other.dataset_name;
    this->dims[0] = other.dims[0];
    this->dims[1] = other.dims[1];
    
    this->iFrame = other.iFrame;
    this->nx = other.nx;
    this->ny = other.ny;
    this->n = other.n;
    
    return *this;
}

PIV::Frame& PIV::Frame::operator=(const PIV::Frame &other)
{
    this->x = std::move(other.x);
    this->y = std::move(other.y);
    this->vx = std::move(other.vx);
    this->vy = std::move(other.vy);
    
    this->dataset_name = other.dataset_name;
    this->dims[0] = other.dims[0];
    this->dims[1] = other.dims[1];
    
    this->iFrame = other.iFrame;
    this->nx = other.nx;
    this->ny = other.ny;
    this->n = other.n;
    
    return *this;
}

PIV::Frame::Frame(const Frame& other)
{
    this->x = std::move(other.x);
    this->y = std::move(other.y);
    this->vx = std::move(other.vx);
    this->vy = std::move(other.vy);
    
    this->dataset_name = other.dataset_name;
    this->dims[0] = other.dims[0];
    this->dims[1] = other.dims[1];
    
    this->iFrame = other.iFrame;
    this->nx = other.nx;
    this->ny = other.ny;
    this->n = other.n;
}

PIV::Frame::~Frame()
{
    this->x.clear();
    this->y.clear();
    this->vx.clear();
    this->vy.clear();
}


void PIV::Frame::unload(H5::H5File& file, std::string groupName)
{
    Group group = file.createGroup(groupName);
    // nx
    DataSpace attr_dataspace = DataSpace(H5S_SCALAR);
    Attribute attr = group.createAttribute("nx", PredType::NATIVE_INT,
                                            attr_dataspace);
    attr.write(PredType::NATIVE_INT, &(this->nx));
    // ny
    attr_dataspace = DataSpace(H5S_SCALAR);
    attr = group.createAttribute("ny", PredType::NATIVE_INT,
                                           attr_dataspace);
    attr.write(PredType::NATIVE_INT, &(this->ny));
    
    hsize_t dims[2] = {static_cast<hsize_t>(this->n), 1};
    // x
    {        
        DataSpace dataspace(1, dims);
        auto dataset = DataSet(file.createDataSet(groupName+"/x",
                                                  PredType::NATIVE_DOUBLE,
                                                  dataspace));
        dataset.write(&(this->x[0]), PredType::NATIVE_DOUBLE);
    }
    
    // y
    {        
        DataSpace dataspace(1, dims);
        auto dataset = DataSet(file.createDataSet(groupName+"/y",
                                                  PredType::NATIVE_DOUBLE,
                                                  dataspace));
        dataset.write(&(this->y[0]), PredType::NATIVE_DOUBLE);
    }
    
    // vx
    {        
        DataSpace dataspace(1, dims);
        auto dataset = DataSet(file.createDataSet(groupName+"/vx",
                                                  PredType::NATIVE_DOUBLE,
                                                  dataspace));
        dataset.write(&(this->vx[0]), PredType::NATIVE_DOUBLE);
    }
    // vy
    {        
        DataSpace dataspace(1, dims);
        auto dataset = DataSet(file.createDataSet(groupName+"/vy",
                                                  PredType::NATIVE_DOUBLE,
                                                  dataspace));
        dataset.write(&(this->vy[0]), PredType::NATIVE_DOUBLE);
    }
}

int PIV::Frame::findWallLocation()
{
    std::vector<double> values(this->ny, 0.0);
    std::vector<double> absV(this->n, 0.0);
    
    for (int i=0; i<this->n; ++i)
    {
        absV[i] = sqrt(this->vx[i]*this->vx[i] +
                       this->vy[i]*this->vy[i]);
    }
    
    for (int iRow=0; iRow<this->ny; ++iRow)
    {
//        for (auto it=absV.begin()+iRow*this->nx;
//             absV.begin()+iRow*this->nx != absV.end();
//             it++)
//        {
//            values[iRow] += *it;
//        }
        values[iRow] = std::accumulate(absV.begin()+iRow*this->nx,
                                       absV.begin()+(iRow+1)*this->nx - 1,
                                       0);
    }
    
    // find minimum
    double minimum = 1e10;
    int imin = 0;
    for (int iRow=0; iRow<this->ny; ++iRow)
    {
        if (values[iRow] < minimum)
        {
            minimum = values[iRow];
            imin = iRow;
        }
    }
    
    return imin;
}













