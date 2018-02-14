#-------------------------------------------------
#
# Project created by QtCreator 2015-08-01T20:41:23
# Jonti Olds
#
#-------------------------------------------------

QT       += multimedia core network gui svg sql

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets  printsupport

TARGET = JAERO
TEMPLATE = app

#message("QT_ARCH is \"$$QT_ARCH\"");
contains(QT_ARCH, i386) {
    #message("32-bit")
    DEFINES += kiss_fft_scalar=float
} else {
    #message("64-bit")
    DEFINES += kiss_fft_scalar=double
}

SOURCES += main.cpp\
        mainwindow.cpp \
    coarsefreqestimate.cpp \
    DSP.cpp \
    mskdemodulator.cpp \
    ../qcustomplot/qcustomplot.cpp \
    audiomskdemodulator.cpp \
    gui_classes/console.cpp \
    gui_classes/qscatterplot.cpp \
    gui_classes/qspectrumdisplay.cpp \
    gui_classes/qled.cpp \
    ../kiss_fft130/kiss_fft.c \
    fftwrapper.cpp \
    fftrwrapper.cpp \
    ../kiss_fft130/kiss_fftr.c \
    gui_classes/textinputwidget.cpp \
    gui_classes/settingsdialog.cpp \
    aerol.cpp \
    gui_classes/planelog.cpp \
    downloadmanager.cpp \
    databasetext.cpp \
    oqpskdemodulator.cpp \
    audiooqpskdemodulator.cpp \
    burstoqpskdemodulator.cpp \
    audioburstoqpskdemodulator.cpp \
    arincparse.cpp \
    tcpserver.cpp \
    sbs1.cpp \
    tcpclient.cpp \
    burstmskdemodulator.cpp \
    audioburstmskdemodulator.cpp \
    jconvolutionalcodec.cpp \
    sdr.cpp \
    ../kiss_fft130/kiss_fastfir_complex.c \
    ../kiss_fft130/kiss_fastfir_real.c


HEADERS  += mainwindow.h \
    coarsefreqestimate.h \
    DSP.h \
    mskdemodulator.h \
    ../qcustomplot/qcustomplot.h \
    audiomskdemodulator.h \
    gui_classes/console.h \
    gui_classes/qscatterplot.h \
    gui_classes/qspectrumdisplay.h \
    gui_classes/qled.h \
    ../kiss_fft130/_kiss_fft_guts.h \
    ../kiss_fft130/kiss_fft.h \
    fftwrapper.h \
    fftrwrapper.h \
    ../kiss_fft130/kiss_fftr.h \
    gui_classes/textinputwidget.h \
    gui_classes/settingsdialog.h \
    aerol.h \
    gui_classes/planelog.h \
    downloadmanager.h \
    databasetext.h \
    oqpskdemodulator.h \
    audiooqpskdemodulator.h \
    burstoqpskdemodulator.h \
    audioburstoqpskdemodulator.h \
    arincparse.h \
    tcpserver.h \
    sbs1.h \
    tcpclient.h \
    burstmskdemodulator.h \
    audioburstmskdemodulator.h \
    jconvolutionalcodec.h \
    sdr.h \
    ../kiss_fft130/kiss_fastfir_complex.h \
    ../kiss_fft130/kiss_fastfir_real.h


FORMS    += mainwindow.ui \
    gui_classes/settingsdialog.ui \
    gui_classes/planelog.ui

RESOURCES += \
    jaero.qrc

DISTFILES += \
    LICENSE \
    ../kiss_fft130/TIPS \
    ../kiss_fft130/CHANGELOG \
    ../kiss_fft130/COPYING \
    ../kiss_fft130/README \
    ../qcustomplot/changelog.txt \
    ../qcustomplot/GPL.txt \
    ../README.md \
    ../images/screenshot-win-main.png \
    ../images/screenshot-win-planelog.png

win32 {
RC_FILE = jaero.rc

}

win32 {
#on windows the libcorrect dlls are here
INCLUDEPATH +=../libcorrect/include
contains(QT_ARCH, i386) {
    #message("32-bit")
    LIBS += -L$$PWD/../libcorrect/bin/32
} else {
    #message("64-bit")
    LIBS += -L$$PWD/../libcorrect/bin/64
}
LIBS += -llibcorrect
}

# remove possible other optimization flags
#QMAKE_CXXFLAGS_RELEASE -= -O
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3

# add the desired -O3 if not present
#QMAKE_CXXFLAGS_RELEASE *= -O3


#sdr includes
win32 {
#on windows this is where I put the headers of the 3rd party libraries that are linked to.
    HEADERS  += ../librtlsdr/include/rtl-sdr.h \
    ../librtlsdr/include/rtl-sdr_export.h
}
unix {
# on the pi Qt creator doesn't seem know about /usr/include and /usr/local/include so I've added them here so it knows
    INCLUDEPATH += /usr/include \
    /usr/local/include
#for me they are saved here but you dont need this
    HEADERS  += /usr/include/rtl-sdr.h \
    /usr/include/rtl-sdr_export.h
}

# here we have librtlsdr libs
win32 {
#on windows the libsdr dlls are here
INCLUDEPATH +=../librtlsdr/include
contains(QT_ARCH, i386) {
    #message("32-bit")
    LIBS += -L$$PWD/../librtlsdr/32
} else {
    #message("64-bit")
    LIBS += -L$$PWD/../librtlsdr/64
}
}
LIBS += -lrtlsdr
