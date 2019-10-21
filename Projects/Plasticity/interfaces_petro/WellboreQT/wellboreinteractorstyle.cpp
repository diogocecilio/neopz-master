#include "wellboreinteractorstyle.h"

WellboreInteractorStyle::WellboreInteractorStyle()
{
}

void WellboreInteractorStyle::OnLeftButtonDown()
{
    vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
}

void WellboreInteractorStyle::OnLeftButtonUp()
{
    vtkInteractorStyleTrackballCamera::OnMiddleButtonUp();
}

vtkStandardNewMacro(WellboreInteractorStyle)
