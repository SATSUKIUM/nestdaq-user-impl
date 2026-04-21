#ifndef INFIXTORP_H
#define INFIXTORP_H

// infixtorp.h //

#include <string>
#include <memory>
struct AstNode {
    std::string value;
    std::shared_ptr<AstNode> left, right;
  public:
    std::string to_rpn() {
        std::string rpn;
        if (left) {
            rpn += left->to_rpn();
        }
        if (right) {
            rpn += right->to_rpn();
        }
        rpn += value + " ";
        return rpn;
    }
};

std::shared_ptr<AstNode> infixtorp(const std::string& input);

#endif // INFIXTORP_H