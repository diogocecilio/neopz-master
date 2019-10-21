#ifndef WELLBOREINTERACTORSTYLE_H
#define WELLBOREINTERACTORSTYLE_H

#include <vtkObjectFactory.h>
#include <vtkInteractorStyleTrackballCamera.h>

class WellboreInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
    vtkTypeMacro(WellboreInteractorStyle,vtkInteractorStyleTrackballCamera)

    static WellboreInteractorStyle* New();

    WellboreInteractorStyle();

    virtual void OnLeftButtonDown();
    virtual void OnLeftButtonUp();
};

#endif // WELLBOREINTERACTORSTYLE_H
