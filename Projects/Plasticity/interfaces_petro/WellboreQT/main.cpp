#include <QObject>
#include <QtGui/QApplication>
#include "mainwindow.h"

#include <QSplashScreen>
#include <QTimer>
#include <QObject>
#include <QTranslator>

#ifdef LOG4CXX
#include "pzlog.h"
#endif

int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif

//MAINWINDOW NO SPLASH CODE
    QApplication a(argc, argv);

    QTranslator translator;
    translator.load("tr_ptbr");
    a.installTranslator(&translator);

    MainWindow w;
    w.showMaximized();

    return a.exec();

//SPLASH WINDOW CODE
//    QApplication a(argc, argv);
//    QPixmap pixmap(":imgs/LogoPetro.png");
//    QSplashScreen *splash = new QSplashScreen(pixmap);
//    splash->showMessage("Loading...");
//    splash->show();

//    MainWindow w;
//    splash->finish(&w);

//    QTimer t;
//    t.setInterval(3000);
//    t.setSingleShot(true);
//    QObject::connect(&t,SIGNAL(timeout()),&w,SLOT(showMaximized()));
//    QObject::connect(&t,SIGNAL(timeout()),splash,SLOT(close()));
//    t.start();

    return a.exec();
}
