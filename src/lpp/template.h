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

#ifndef LOCI_TEMPLATE_H
#define LOCI_TEMPLATE_H

#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Loci {

class TemplateContext {
public:
  TemplateContext();
  TemplateContext(const TemplateContext& other);
  TemplateContext(TemplateContext&& other) noexcept;
  TemplateContext& operator=(const TemplateContext& other);
  TemplateContext& operator=(TemplateContext&& other) noexcept;
  ~TemplateContext();

  void set(std::string key, std::string value);
  void set(std::string key, const char* value);
  void set(std::string key, bool value);
  void set(std::string key, char value);
  void set(std::string key, short value);
  void set(std::string key, int value);
  void set(std::string key, long value);
  void set(std::string key, long long value);
  void set(std::string key, unsigned short value);
  void set(std::string key, unsigned int value);
  void set(std::string key, unsigned long value);
  void set(std::string key, unsigned long long value);
  void set(std::string key, float value);
  void set(std::string key, double value);
  void set(std::string key, long double value);

  void set_object(std::string key, TemplateContext object);
  void set_array(std::string key, std::vector<TemplateContext> items);

  std::optional<std::string> resolve(std::string_view path) const;
  const std::vector<TemplateContext>* array_at(std::string_view path) const;
  bool has_path(std::string_view path) const;
  bool is_true(std::string_view path) const;

  static TemplateContext from_flat(
    const std::unordered_map<std::string, std::string>& values);

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;

  template<typename T>
  static bool value_is_true(const T& value);

  bool is_empty() const;
};

struct RenderOptions {
  bool keep_missing_placeholders = false;
};

class TemplateEngine {
public:
  explicit TemplateEngine(RenderOptions options = {});

  std::string render(std::string_view template_text, const TemplateContext& context) const;

  std::string render_file(const std::string& path, const TemplateContext& context) const;

private:
  RenderOptions options_;
};

}  // namespace Loci

#endif  // LOCI_TEMPLATE_H
