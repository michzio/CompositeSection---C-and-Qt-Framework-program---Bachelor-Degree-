#include "mouseeventshandler.h"

MouseEventsHandler::MouseEventsHandler(QObject *parent) :
    QObject(parent)
{
}


bool MouseEventsHandler::eventFilter(QObject *obj, QEvent *event)
 {
     if (event->type() == QEvent::MouseMove) {
         fprintf(stderr, "Ruch myszki\n");
         return true;
     } else if(event->type() == QEvent::MouseButtonPress) {
         fprintf(stderr, "Wcisnieto przycisk myszy\n");
         return true;
     } else if(event->type() == QEvent::MouseButtonRelease) {
         fprintf(stderr, "Puszczono przycisk myszy\n");
         return true;
     } else {
         // standard event processing
         return QObject::eventFilter(obj, event);
     }
 }
