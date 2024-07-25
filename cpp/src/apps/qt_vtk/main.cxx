#include "MainTigerWindow.h"
#include <QApplication>
#include <QSurfaceFormat>
#include <QVTKOpenGLNativeWidget.h>

int main( int argc, char** argv )
{
	// needed to ensure appropriate OpenGL context is created for VTK rendering
	QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
	
	// QT Stuff
	QApplication app( argc, argv );
	
	// Start the main TIGER window and show it
	MainTigerWindow MainTigerWindow;
	MainTigerWindow.show();
	
	// 
	return app.exec();
}
