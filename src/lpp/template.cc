//#############################################################################
//#
//# Copyright 2008-2025, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

#include "template.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>

namespace Loci {

namespace {

std::string head_segment(std::string_view path) {
  const auto dot = path.find('.');
  if (dot == std::string_view::npos) {
    return std::string(path);
  }
  return std::string(path.substr(0, dot));
}

std::string_view tail_segment(std::string_view path) {
  const auto dot = path.find('.');
  if (dot == std::string_view::npos) {
    return {};
  }
  return path.substr(dot + 1);
}

std::string lower_copy(std::string_view value) {
  std::string out(value);
  std::transform(out.begin(), out.end(), out.begin(),
         [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return out;
}

bool is_true_string(std::string_view value) {
  if (value.empty()) {
    return false;
  }

  const std::string lower = lower_copy(value);
  if (lower == "false" || lower == "0" || lower == "no" || lower == "off") {
    return false;
  }
  return true;
}

std::string trim_copy(const std::string& s) {
  const auto start = s.find_first_not_of(" \t\r\n");
  if (start == std::string::npos) {
    return {};
  }
  const auto end = s.find_last_not_of(" \t\r\n");
  return s.substr(start, end - start + 1);
}

struct PathRef {
  bool from_root = false;
  std::string_view path;
};

std::string_view trim_leading_path(std::string_view path) {
  const auto start = path.find_first_not_of(" \t");
  if (start == std::string_view::npos) {
    return {};
  }
  return path.substr(start);
}

PathRef parse_path_ref(std::string_view path) {
  if (path.size() >= 2 && path[0] == ':' && path[1] == ':') {
    return {true, trim_leading_path(path.substr(2))};
  }
  return {false, path};
}

bool is_path_char(char c) {
  return std::isalnum(static_cast<unsigned char>(c)) != 0 || c == '_' || c == '.' || c == '@';
}

[[noreturn]] void throw_condition_error(const std::string& message) {
  throw std::runtime_error("Invalid condition expression: " + message);
}

struct RenderState;

enum class CondTokenKind { Path, And, Or, Not, LParen, RParen, End };

struct CondLexer {
  std::string_view input;
  std::size_t pos = 0;
  CondTokenKind kind = CondTokenKind::End;
  std::string path;

  explicit CondLexer(std::string_view input) : input(input) {
    advance();
  }

  void skip_space() {
    while (pos < input.size() &&
       (input[pos] == ' ' || input[pos] == '\t')) {
      ++pos;
    }
  }

  void advance() {
    skip_space();
    if (pos >= input.size()) {
      kind = CondTokenKind::End;
      path.clear();
      return;
    }

    const char c = input[pos];
    if (c == '(') {
      ++pos;
      kind = CondTokenKind::LParen;
      path.clear();
      return;
    }
    if (c == ')') {
      ++pos;
      kind = CondTokenKind::RParen;
      path.clear();
      return;
    }
    if (c == '!') {
      ++pos;
      kind = CondTokenKind::Not;
      path.clear();
      return;
    }
    if (c == '&' && pos + 1 < input.size() && input[pos + 1] == '&') {
      pos += 2;
      kind = CondTokenKind::And;
      path.clear();
      return;
    }
    if (c == '|' && pos + 1 < input.size() && input[pos + 1] == '|') {
      pos += 2;
      kind = CondTokenKind::Or;
      path.clear();
      return;
    }
    if (c == '&' || c == '|') {
      throw_condition_error("expected '&&' or '||' operator");
    }

    const std::size_t start = pos;
    if (pos + 1 < input.size() && input[pos] == ':' && input[pos + 1] == ':') {
      pos += 2;
      skip_space();
      if (pos >= input.size() || !is_path_char(input[pos])) {
        throw_condition_error("expected path after '::'");
      }
    }
    if (pos >= input.size() || !is_path_char(input[pos])) {
      throw_condition_error(std::string("expected path, found '") + c + "'");
    }
    while (pos < input.size() && is_path_char(input[pos])) {
      ++pos;
    }
    if (pos == start) {
      throw_condition_error("expected path");
    }

    kind = CondTokenKind::Path;
    path = std::string(input.substr(start, pos - start));
  }

  void expect(CondTokenKind expected) {
    if (kind != expected) {
      throw_condition_error("unexpected token in condition");
    }
    advance();
  }
};

struct ConditionEvaluator;

struct TemplateNode {
  enum class Type { Text, Variable, IfBlock, EachBlock };

  struct ElifBranch {
    std::string condition;
    std::vector<TemplateNode> body;
  };

  Type type = Type::Text;
  std::string value;
  std::vector<TemplateNode> then_branch;
  std::vector<ElifBranch> elif_branches;
  std::vector<TemplateNode> else_branch;
  std::vector<TemplateNode> body;
};

using NodeList = std::vector<TemplateNode>;

enum class TokenKind {
  Text,
  Variable,
  IfOpen,
  ElifOpen,
  Else,
  IfClose,
  EachOpen,
  EachClose,
  End
};

struct SourceLocation {
  std::size_t line = 1;
  std::size_t column = 1;
};

struct Token {
  TokenKind kind = TokenKind::End;
  std::string value;
  SourceLocation location;
};

std::string token_kind_label(TokenKind kind) {
  switch (kind) {
    case TokenKind::Text:
      return "text";
    case TokenKind::Variable:
      return "variable";
    case TokenKind::IfOpen:
      return "#if";
    case TokenKind::ElifOpen:
      return "#else if";
    case TokenKind::Else:
      return "#else";
    case TokenKind::IfClose:
      return "/if";
    case TokenKind::EachOpen:
      return "#each";
    case TokenKind::EachClose:
      return "/each";
    case TokenKind::End:
      return "end of template";
  }
  return "unknown token";
}

[[noreturn]] void throw_parse_error(const SourceLocation& location, const std::string& message) {
  std::ostringstream out;
  out << message << " at line " << location.line << ", column " << location.column;
  throw std::runtime_error(out.str());
}

[[noreturn]] void throw_parse_error(const Token& token, const std::string& message) {
  std::ostringstream out;
  out << message << " at line " << token.location.line << ", column " << token.location.column;
  if (token.kind != TokenKind::End && token.kind != TokenKind::Text) {
    out << " (found " << token_kind_label(token.kind);
    if (!token.value.empty()) {
      out << " '" << token.value << "'";
    }
    out << ")";
  }
  throw std::runtime_error(out.str());
}

[[noreturn]] void throw_unknown_control(const std::string& content, SourceLocation location) {
  std::ostringstream out;
  out << "Unrecognized control statement '" << content << "'";
  if (content == "/fi") {
    out << " (did you mean '/if'?)";
  }
  throw_parse_error(location, out.str());
}

class Tokenizer {
public:
  explicit Tokenizer(std::string_view input) : input_(input) {}

  Token next() {
    return lex_next();
  }

private:
  std::string_view input_;
  std::size_t pos_ = 0;
  std::size_t line_ = 1;
  std::size_t column_ = 1;

  SourceLocation current_location() const {
    return {line_, column_};
  }

  void advance_to(std::size_t new_pos) {
    while (pos_ < new_pos && pos_ < input_.size()) {
      if (input_[pos_] == '\n') {
        ++line_;
        column_ = 1;
      } else {
        ++column_;
      }
      ++pos_;
    }
  }

  Token lex_next() {
    if (pos_ >= input_.size()) {
      return {TokenKind::End, {}, current_location()};
    }

    const std::size_t open = input_.find("$$", pos_);
    if (open == std::string_view::npos) {
      const SourceLocation location = current_location();
      Token token{TokenKind::Text, std::string(input_.substr(pos_)), location};
      advance_to(input_.size());
      return token;
    }

    if (open > pos_) {
      const SourceLocation location = current_location();
      Token token{TokenKind::Text, std::string(input_.substr(pos_, open - pos_)), location};
      advance_to(open);
      return token;
    }

    const SourceLocation location = current_location();
    const std::size_t content_start = open + 2;
    const std::size_t close = input_.find("$$", content_start);
    if (close == std::string_view::npos) {
      Token token{TokenKind::Text, std::string(input_.substr(pos_)), location};
      advance_to(input_.size());
      return token;
    }

    const std::string content =
      trim_copy(std::string(input_.substr(content_start, close - content_start)));
    advance_to(close + 2);
    return classify(content, location);
  }

  static Token classify(const std::string& content, SourceLocation location) {
    if (content == "/if") {
      return {TokenKind::IfClose, {}, location};
    }
    if (content == "/each") {
      return {TokenKind::EachClose, {}, location};
    }
    if (content == "#else") {
      return {TokenKind::Else, {}, location};
    }
    if (content.rfind("#else if ", 0) == 0) {
      const std::string condition = trim_copy(content.substr(9));
      if (condition.empty()) {
        throw_parse_error(location, "Expected condition after '#else if'");
      }
      return {TokenKind::ElifOpen, condition, location};
    }
    if (content.rfind("#elseif ", 0) == 0) {
      const std::string condition = trim_copy(content.substr(8));
      if (condition.empty()) {
        throw_parse_error(location, "Expected condition after '#elseif'");
      }
      return {TokenKind::ElifOpen, condition, location};
    }
    if (content.rfind("#if ", 0) == 0) {
      const std::string condition = trim_copy(content.substr(4));
      if (condition.empty()) {
        throw_parse_error(location, "Expected condition after '#if'");
      }
      return {TokenKind::IfOpen, condition, location};
    }
    if (content == "#if") {
      throw_parse_error(location, "Expected condition after '#if'");
    }
    if (content.rfind("#each ", 0) == 0) {
      const std::string list = trim_copy(content.substr(6));
      if (list.empty()) {
        throw_parse_error(location, "Expected list path after '#each'");
      }
      return {TokenKind::EachOpen, list, location};
    }
    if (content == "#each") {
      throw_parse_error(location, "Expected list path after '#each'");
    }

    if (!content.empty() && (content[0] == '#' || content[0] == '/')) {
      throw_unknown_control(content, location);
    }

    return {TokenKind::Variable, content, location};
  }
};

class Parser {
public:
  explicit Parser(std::string_view input) : tokenizer_(input) {
    current_ = tokenizer_.next();
  }

  NodeList parse_document() {
    NodeList nodes =
      parse_segments({TokenKind::IfClose, TokenKind::EachClose, TokenKind::End});
    if (current_.kind != TokenKind::End) {
      throw_parse_error(current_, "Unexpected token after document end");
    }
    return nodes;
  }

private:
  Tokenizer tokenizer_;
  Token current_;

  void advance() {
    current_ = tokenizer_.next();
  }

  bool at(TokenKind kind) const {
    return current_.kind == kind;
  }

  bool at_end() const {
    return current_.kind == TokenKind::End;
  }

  void expect(TokenKind kind, const char* message) {
    if (!at(kind)) {
      throw_parse_error(current_, message);
    }
    advance();
  }

  NodeList parse_segments(std::initializer_list<TokenKind> stop_tokens) {
    NodeList nodes;
    while (!at_end()) {
      bool stop = false;
      for (TokenKind kind : stop_tokens) {
        if (at(kind)) {
          stop = true;
          break;
        }
      }
      if (stop) {
        break;
      }
      nodes.push_back(parse_segment());
    }
    return nodes;
  }

  TemplateNode parse_segment() {
    switch (current_.kind) {
      case TokenKind::Text: {
        TemplateNode node;
        node.type = TemplateNode::Type::Text;
        node.value = std::move(current_.value);
        advance();
        return node;
      }
      case TokenKind::Variable: {
        TemplateNode node;
        node.type = TemplateNode::Type::Variable;
        node.value = std::move(current_.value);
        advance();
        return node;
      }
      case TokenKind::IfOpen:
        return parse_if_block();
      case TokenKind::EachOpen:
        return parse_each_block();
      default:
        throw_parse_error(current_, "Unexpected token in template segment");
    }
  }

  TemplateNode parse_if_block() {
    TemplateNode node;
    node.type = TemplateNode::Type::IfBlock;
    node.value = std::move(current_.value);
    advance();

    const std::initializer_list<TokenKind> branch_stops = {
      TokenKind::IfClose, TokenKind::Else, TokenKind::ElifOpen};

    node.then_branch = parse_segments(branch_stops);
    while (at(TokenKind::ElifOpen)) {
      TemplateNode::ElifBranch branch;
      branch.condition = std::move(current_.value);
      advance();
      branch.body = parse_segments(branch_stops);
      node.elif_branches.push_back(std::move(branch));
    }
    if (at(TokenKind::Else)) {
      advance();
      node.else_branch = parse_segments({TokenKind::IfClose});
    }
    expect(TokenKind::IfClose, "Expected $$/if$$");
    return node;
  }

  TemplateNode parse_each_block() {
    TemplateNode node;
    node.type = TemplateNode::Type::EachBlock;
    node.value = std::move(current_.value);
    advance();

    node.body = parse_segments({TokenKind::EachClose});
    expect(TokenKind::EachClose, "Expected $$/each$$");
    return node;
  }
};

struct RenderState {
  const TemplateContext& root;
  const RenderOptions& options;
  std::vector<const TemplateContext*> scopes;
  std::string result;

  const TemplateContext* owning_context(std::string_view path) const {
    for (auto it = scopes.rbegin(); it != scopes.rend(); ++it) {
      if ((*it)->has_path(path)) {
        return *it;
      }
    }
    if (root.has_path(path)) {
      return &root;
    }
    return nullptr;
  }

  std::optional<std::string> lookup(std::string_view raw_path) const {
    const PathRef ref = parse_path_ref(raw_path);
    if (ref.path.empty()) {
      return std::nullopt;
    }
    if (ref.from_root) {
      return root.resolve(ref.path);
    }
    if (const TemplateContext* ctx = owning_context(ref.path)) {
      return ctx->resolve(ref.path);
    }
    return std::nullopt;
  }

  bool path_is_true(std::string_view raw_path) const {
    const PathRef ref = parse_path_ref(raw_path);
    if (ref.path.empty()) {
      return false;
    }
    if (ref.from_root) {
      return root.is_true(ref.path);
    }
    if (const TemplateContext* ctx = owning_context(ref.path)) {
      return ctx->is_true(ref.path);
    }
    return false;
  }

  bool is_true(std::string_view condition) const;

  const std::vector<TemplateContext>* array_at(std::string_view raw_path) const {
    const PathRef ref = parse_path_ref(raw_path);
    if (ref.path.empty()) {
      return nullptr;
    }
    if (ref.from_root) {
      return root.array_at(ref.path);
    }
    if (const TemplateContext* ctx = owning_context(ref.path)) {
      return ctx->array_at(ref.path);
    }
    return nullptr;
  }

  void append_variable(std::string_view path) {
    if (const auto value = lookup(path)) {
      result += *value;
      return;
    }

    if (options.keep_missing_placeholders) {
      result += "$$";
      result += path;
      result += "$$";
    }
  }
};

struct ConditionEvaluator {
  const RenderState& state;
  CondLexer lex;

  ConditionEvaluator(const RenderState& state, std::string_view input)
    : state(state), lex(input) {}

  bool parse_or();
  bool parse_and();
  bool parse_unary();
  bool parse_primary();
  bool evaluate();
};

bool ConditionEvaluator::parse_or() {
  bool value = parse_and();
  while (lex.kind == CondTokenKind::Or) {
    lex.advance();
    const bool rhs = parse_and();
    value = value || rhs;
  }
  return value;
}

bool ConditionEvaluator::parse_and() {
  bool value = parse_unary();
  while (lex.kind == CondTokenKind::And) {
    lex.advance();
    const bool rhs = parse_unary();
    value = value && rhs;
  }
  return value;
}

bool ConditionEvaluator::parse_unary() {
  if (lex.kind == CondTokenKind::Not) {
    lex.advance();
    return !parse_unary();
  }
  return parse_primary();
}

bool ConditionEvaluator::parse_primary() {
  if (lex.kind == CondTokenKind::Path) {
    const std::string path = std::move(lex.path);
    lex.advance();
    return state.path_is_true(path);
  }
  if (lex.kind == CondTokenKind::LParen) {
    lex.advance();
    const bool value = parse_or();
    lex.expect(CondTokenKind::RParen);
    return value;
  }
  throw_condition_error("expected path or '('");
}

bool ConditionEvaluator::evaluate() {
  if (lex.kind == CondTokenKind::End) {
    throw_condition_error("empty condition");
  }
  const bool value = parse_or();
  if (lex.kind != CondTokenKind::End) {
    throw_condition_error("unexpected token after expression");
  }
  return value;
}

bool RenderState::is_true(std::string_view condition) const {
  return ConditionEvaluator(*this, condition).evaluate();
}

void render_nodes(const NodeList& nodes, RenderState& state);

void render_node(const TemplateNode& node, RenderState& state) {
  switch (node.type) {
    case TemplateNode::Type::Text:
      state.result += node.value;
      break;

    case TemplateNode::Type::Variable:
      state.append_variable(node.value);
      break;

    case TemplateNode::Type::IfBlock:
      if (state.is_true(node.value)) {
        render_nodes(node.then_branch, state);
      } else {
        bool matched = false;
        for (const auto& elif : node.elif_branches) {
          if (state.is_true(elif.condition)) {
            render_nodes(elif.body, state);
            matched = true;
            break;
          }
        }
        if (!matched) {
          render_nodes(node.else_branch, state);
        }
      }
      break;

    case TemplateNode::Type::EachBlock: {
      const auto* items = state.array_at(node.value);
      if (items == nullptr) {
        break;
      }
      const int size = static_cast<int>(items->size());
      for (std::size_t i = 0; i < items->size(); ++i) {
        TemplateContext loop_ctx = (*items)[i];
        loop_ctx.set("@index", static_cast<int>(i));
        loop_ctx.set("@size", size);
        state.scopes.push_back(&loop_ctx);
        render_nodes(node.body, state);
        state.scopes.pop_back();
      }
      break;
    }
  }
}

void render_nodes(const NodeList& nodes, RenderState& state) {
  for (const auto& node : nodes) {
    render_node(node, state);
  }
}

std::string parse_and_render(std::string_view template_text,
              const TemplateContext& context,
              const RenderOptions& options) {
  Parser parser(template_text);
  const NodeList document = parser.parse_document();

  RenderState state{context, options, {}, {}};
  render_nodes(document, state);
  return state.result;
}

using ContextValue = std::variant<
  std::string,
  bool,
  char,
  short,
  int,
  long,
  long long,
  unsigned short,
  unsigned int,
  unsigned long,
  unsigned long long,
  float,
  double,
  long double,
  TemplateContext,
  std::vector<TemplateContext>>;

std::string value_as_string(const ContextValue& value) {
  return std::visit(
    [](const auto& v) -> std::string {
      using T = std::decay_t<decltype(v)>;
      if constexpr (std::is_same_v<T, std::string>) {
        return v;
      } else if constexpr (std::is_same_v<T, bool>) {
        return v ? "1" : "0";
      } else if constexpr (std::is_same_v<T, char>) {
        return std::string(1, v);
      } else if constexpr (std::is_same_v<T, TemplateContext> ||
                 std::is_same_v<T, std::vector<TemplateContext>>) {
        return {};
      } else {
        return std::to_string(v);
      }
    },
    value);
}

}  // namespace

template<>
bool TemplateContext::value_is_true<ContextValue>(const ContextValue& value) {
  return std::visit(
    [](const auto& v) -> bool {
      using T = std::decay_t<decltype(v)>;
      if constexpr (std::is_same_v<T, std::string>) {
        return is_true_string(v);
      } else if constexpr (std::is_same_v<T, bool>) {
        return v;
      } else if constexpr (std::is_same_v<T, char>) {
        return v != '\0';
      } else if constexpr (std::is_same_v<T, std::vector<TemplateContext>>) {
        return !v.empty();
      } else if constexpr (std::is_same_v<T, TemplateContext>) {
        return !v.is_empty();
      } else {
        return v != 0;
      }
    },
    value);
}

struct TemplateContext::Impl {
  std::unordered_map<std::string, ContextValue> members_;
};

TemplateContext::TemplateContext()
  : impl_(std::make_unique<Impl>()) {}

TemplateContext::TemplateContext(const TemplateContext& other)
  : impl_(other.impl_ ? std::make_unique<Impl>(*other.impl_)
                      : std::make_unique<Impl>()) {}

TemplateContext::TemplateContext(TemplateContext&& other) noexcept = default;

TemplateContext& TemplateContext::operator=(const TemplateContext& other) {
  if (this != &other) {
    impl_ = other.impl_ ? std::make_unique<Impl>(*other.impl_)
                        : std::make_unique<Impl>();
  }
  return *this;
}

TemplateContext& TemplateContext::operator=(TemplateContext&& other) noexcept = default;

TemplateContext::~TemplateContext() = default;

void TemplateContext::set(std::string key, std::string value) {
  impl_->members_[std::move(key)] = std::move(value);
}

void TemplateContext::set(std::string key, const char* value) {
  set(std::move(key), std::string(value));
}

void TemplateContext::set(std::string key, bool value) {
  impl_->members_[std::move(key)] = value;
}

void TemplateContext::set(std::string key, char value) {
  impl_->members_[std::move(key)] = value;
}

void TemplateContext::set(std::string key, short value) {
  impl_->members_[std::move(key)] = value;
}

void TemplateContext::set(std::string key, int value) {
  impl_->members_[std::move(key)] = value;
}

void TemplateContext::set(std::string key, long value) {
  impl_->members_[std::move(key)] = value;
}

void TemplateContext::set(std::string key, long long value) {
  impl_->members_[std::move(key)] = value;
}

void TemplateContext::set(std::string key, unsigned short value) {
  impl_->members_[std::move(key)] = value;
}

void TemplateContext::set(std::string key, unsigned int value) {
  impl_->members_[std::move(key)] = value;
}

void TemplateContext::set(std::string key, unsigned long value) {
  impl_->members_[std::move(key)] = value;
}

void TemplateContext::set(std::string key, unsigned long long value) {
  impl_->members_[std::move(key)] = value;
}

void TemplateContext::set(std::string key, float value) {
  impl_->members_[std::move(key)] = value;
}

void TemplateContext::set(std::string key, double value) {
  impl_->members_[std::move(key)] = value;
}

void TemplateContext::set(std::string key, long double value) {
  impl_->members_[std::move(key)] = value;
}

void TemplateContext::set_object(std::string key, TemplateContext object) {
  impl_->members_[std::move(key)] = std::move(object);
}

void TemplateContext::set_array(std::string key, std::vector<TemplateContext> items) {
  impl_->members_[std::move(key)] = std::move(items);
}

TemplateContext TemplateContext::from_flat(
  const std::unordered_map<std::string, std::string>& values) {
  TemplateContext ctx;
  for (const auto& [key, value] : values) {
    ctx.set(key, value);
  }
  return ctx;
}

std::optional<std::string> TemplateContext::resolve(std::string_view path) const {
  if (path.empty()) {
    return std::nullopt;
  }

  const std::string head = head_segment(path);
  const std::string_view tail = tail_segment(path);

  const auto it = impl_->members_.find(head);
  if (it == impl_->members_.end()) {
    return std::nullopt;
  }

  if (!tail.empty()) {
    if (const auto* object = std::get_if<TemplateContext>(&it->second)) {
      return object->resolve(tail);
    }
    return std::nullopt;
  }

  if (std::holds_alternative<TemplateContext>(it->second) ||
    std::holds_alternative<std::vector<TemplateContext>>(it->second)) {
    return std::nullopt;
  }

  return value_as_string(it->second);
}

const std::vector<TemplateContext>* TemplateContext::array_at(std::string_view path) const {
  if (path.empty()) {
    return nullptr;
  }

  const std::string head = head_segment(path);
  const std::string_view tail = tail_segment(path);

  const auto it = impl_->members_.find(head);
  if (it == impl_->members_.end()) {
    return nullptr;
  }

  if (!tail.empty()) {
    if (const auto* object = std::get_if<TemplateContext>(&it->second)) {
      return object->array_at(tail);
    }
    return nullptr;
  }

  if (const auto* array = std::get_if<std::vector<TemplateContext>>(&it->second)) {
    return array;
  }
  return nullptr;
}

bool TemplateContext::has_path(std::string_view path) const {
  if (path.empty()) {
    return false;
  }

  const std::string head = head_segment(path);
  const std::string_view tail = tail_segment(path);

  const auto it = impl_->members_.find(head);
  if (it == impl_->members_.end()) {
    return false;
  }

  if (tail.empty()) {
    return true;
  }

  if (const auto* object = std::get_if<TemplateContext>(&it->second)) {
    return object->has_path(tail);
  }
  return false;
}

bool TemplateContext::is_empty() const {
  return impl_->members_.empty();
}

bool TemplateContext::is_true(std::string_view path) const {
  if (path.empty()) {
    return false;
  }

  const std::string head = head_segment(path);
  const std::string_view tail = tail_segment(path);

  const auto it = impl_->members_.find(head);
  if (it == impl_->members_.end()) {
    return false;
  }

  if (tail.empty()) {
    return value_is_true<ContextValue>(it->second);
  }

  if (const auto* object = std::get_if<TemplateContext>(&it->second)) {
    return object->is_true(tail);
  }
  return false;
}

TemplateEngine::TemplateEngine(RenderOptions options)
  : options_(std::move(options)) {}

std::string TemplateEngine::render(std::string_view template_text,
                 const TemplateContext& context) const {
  return parse_and_render(template_text, context, options_);
}

std::string TemplateEngine::render_file(const std::string& path,
                    const TemplateContext& context) const {
  std::ifstream file(path);
  if (!file) {
    throw std::runtime_error("Failed to open template file: " + path);
  }

  std::ostringstream buffer;
  buffer << file.rdbuf();
  return render(buffer.str(), context);
}

}  // namespace Loci
