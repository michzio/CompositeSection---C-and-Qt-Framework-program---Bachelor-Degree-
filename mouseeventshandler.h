#ifndef MOUSEEVENTSHANDLER_H
#define MOUSEEVENTSHANDLER_H

#include <QEvent>
#include <QObject>

class MouseEventsHandler : public QObject
{
    Q_OBJECT
public:
    explicit MouseEventsHandler(QObject *parent = 0);
    
protected:
     bool eventFilter(QObject *obj, QEvent *event);
    
};

#endif // MOUSEEVENTSHANDLER_H
