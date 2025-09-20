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
#include <Tools/variable.h>

using std::vector ;
using std::string ;
using std::make_pair ;

using std::cerr ;
using std::endl ;
using std::ostream ;
using std::istream ;

namespace Loci {
  time_ident::time_hierarchy *time_ident::thp = 0 ;
  
  int time_ident::time_hierarchy::add_level(const string &lname,int level) {
    vector<int>::const_iterator vi ;
    for(vi=time_db[level].children.begin();
        vi!=time_db[level].children.end();++vi) {
      if(time_db[*vi].level_name == lname)
	break ;
    }
    if(vi==time_db[level].children.end()) {
      time_db.push_back(time_info(lname,level)) ;
      int n = time_db.size() - 1 ;
      time_db[level].children.push_back(n) ;
      return n ;
    }
    else
      return *vi ;
  } 
   
  time_ident::time_ident(const exprP &exp) {
    create_thp() ;
    id = 0 ;
    if(exp->op == OP_NAME) {
      id = thp->add_level(exp->name,id) ;
      return ;
    }
    if(exp->op != OP_NAME_BRACE) {
      cerr << "syntax error interpreting time list in expression "
	   << exp << endl ;
      return ;
    }
    exprList::const_iterator bgn,end,ii ;
    bgn = exp->expr_list.begin() ;
    end = exp->expr_list.end() ;
    --end ;
    for(ii=bgn;ii!=end;++ii) {
      if(OP_NAME != (*ii)->op) {
	cerr << "syntax error interpreting time list in expression "
	     << exp << endl ;
	cerr << "expecting OP_NAME, got " << *ii << endl ;
	break ;
      }
      else 
	id = thp->add_level((*ii)->name,id) ;
    }
    
    if(OP_NAME == (*end)->op) {
      id = thp->add_level((*ii)->name,id) ;
      return ;
    }
    
    if(OP_PLUS == (*end)->op ||
       OP_MINUS == (*end)->op ||
       OP_ASSIGN == (*end)->op) {
        if((*end)->expr_list.size() != 2) {
            cerr << "syntax error in time label of name expression "
                 << exp << endl ;
            cerr << "error found while parsing " << *end << endl ;
            id = 0 ;
            return ;
        }
        exprP tname = (*end)->expr_list.front() ;
        if(OP_NAME != tname->op) {
            cerr << "syntax error in time label of name expression "
                 << exp << endl ;
            cerr << "time level should be of form NAME <op> int" << endl;
            id = 0 ;
            return ;
        } else {

          id = thp->add_level(tname->name,id) ;
	  return ;
        }
    }
    
    //cerr << "syntax error in time label of name expression "
    //   << exp << endl ;
    id = 0 ;
  }

  time_ident::time_ident(const std::string &lname, const time_ident &t) {
    std::list<std::string> names ;
    time_ident stationary ;
    time_ident current = t ;
    while(current!=stationary) {
      names.push_front(current.level_name()) ;
      current = current.parent() ;
    }
    names.push_front(lname) ;
    std::list<std::string>::const_iterator li ;
    create_thp() ;
    id = 0 ;
    for(li=names.begin();li!=names.end();++li)
      id = thp->add_level(*li,id) ;
  }
  time_ident::time_ident(const time_ident &t, const std::string &lname) {
    create_thp() ;
    id = thp->add_level(lname,t.ident()) ;
  }                                                                    
                                                                    
  
ostream &time_ident::Print(ostream &s) const {
    if(id==0)
      return s ;
    vector<int> v ;
    int l = id ;
    while(l != 0) {
        v.push_back(l) ;
        l = thp->parent(l) ;
    }
    s << thp->name(v.back()) ;
    v.pop_back() ;
    while(v.size() != 0) {
        s << "," << thp->name(v.back()) ;
        v.pop_back() ;
    }
    return s ;
}

bool time_ident::before(const time_ident &t) const {
  int i = t.id ;
  while(i != 0) {
    i = thp->parent(i) ;
    if(i == id)
      return true ;
  }
  return false ;
}

vector<time_ident> time_ident::children() {
  const vector<int> &vc = thp->get_children(id) ;
  vector<time_ident> v ;
  for(size_t i=0;i<vc.size();++i)
    v.push_back(time_ident(vc[i])) ;
  return v ;
}


variable::variable_db *variable::vdb = 0 ;

bool variable::info::operator<(const info &v) const {

    if(tvar != v.tvar)
      return tvar < v.tvar ;
    else if(name != v.name)
      return name < v.name ;
    else if(time_id != v.time_id)
      return time_id < v.time_id ;
    else if(time_id != time_ident() && assign != v.assign)
      return assign < v.assign ;
    else if(time_id != time_ident() && offset != v.offset)
      return offset < v.offset ;
    else if(priority != v.priority)
      return priority < v.priority ;
    else if(namespac != v.namespac)
      return namespac < v.namespac ;
    else
      return v_ids < v.v_ids ;
}

  
  bool variable::info::operator==(const info &v) const {
    return
      tvar    == v.tvar     &&
      assign  == v.assign   &&
      namespac== v.namespac &&
      name    == v.name     &&
      time_id == v.time_id  &&
      offset  == v.offset   &&
      //rename  == v.rename   &&
      v_ids   == v.v_ids ;
  }
  
  variable::variable(const exprP &p) {
    create_vdb() ;
    info v ;
    exprP e = p ;
    if(OP_AT == e->op) {
      exprList l = collect_associative_op(e,OP_AT);
      while(l.begin() != l.end()) {
	exprP s = l.front();
	l.pop_front();
	if(l.begin()==l.end()) {
	  e = s;
	} else {
	  std::string tmp_name ;
	  switch(s->op) {
	  case OP_NAME:
	    //v.rename = true ;
	    v.namespac.push_back(s->name);
	    break ;
	  case OP_SCOPE:
            {
              exprList l = collect_associative_op(s,OP_SCOPE) ;
              
              while(l.begin() != l.end()) {
                exprP s = l.front() ;
                l.pop_front() ;
                if(l.begin() == l.end()) 
                  v.namespac.push_back(s->name) ; //e = s ;
                else 
                  if(OP_NAME == s->op) {
                    v.priority.push_back(s->name) ;
                  } else {
                    cerr << "unable to interpret priority list in expression "
                         << s << endl
                         << "error occured while parsing " << e << endl ;
                  }
              }
            }
	    break ;
	  default:
	    cerr << "unable to interpret namespace list in expression " << s << endl << "error occured while parsing " << e << endl;
	    break;
	  } 
	} 
      } 
    } 
    
    
    if(OP_SCOPE == e->op) {
      exprList l = collect_associative_op(e,OP_SCOPE) ;
      while(l.begin() != l.end()) {
	exprP s = l.front() ;
	l.pop_front() ;
	if(l.begin() == l.end()) 
	  e = s ;
	else 
	  switch(s->op) {
	  case OP_NAME:
	    v.priority.push_back(s->name) ;
	    break ;
	  default:
	    cerr << "unable to interpret priority list in expression "
		 << s << endl
		 << "error occured while parsing " << e << endl ;
	    
	  }
      }
    }
    
    if(e->op == OP_DOLLAR) {
      v.tvar = true ;
      e = e->expr_list.front() ;
    }
    exprP end ;
    switch(e->op) {
    case OP_NAME:
      v.name = e->name ;
      break ;
    case OP_NAME_BRACE:
      v.name = e->name ;
      v.time_id = time_ident(e) ;
      end = e->expr_list.back() ;
      if(OP_PLUS == end->op  ||
	 OP_MINUS == end->op ||
	 OP_ASSIGN == end->op) {
	if(end->expr_list.size() != 2) {
	  cerr << "syntax error in time label of name expression "
	       << p << endl ;
	  cerr << "error found while parsing " << end << endl ;
	  break ;
	}
	exprP tname = end->expr_list.front() ;
	int time_offset = 0 ;
	if(OP_NAME != tname->op) {
	  cerr << "syntax error in time label of name expression "
	       << p << endl ;
	  cerr << "time level should be of form NAME <op> int" << endl;
	  break ;
	} else
	  if(OP_INT == end->expr_list.back()->op) {
	    time_offset = end->expr_list.back()->int_val ;
	  } else if((OP_UNARY_PLUS  == end->expr_list.back()->op ||
		     OP_UNARY_MINUS == end->expr_list.back()->op) &&
		    OP_INT == end->expr_list.back()->expr_list.front()->op){
	    time_offset = end->expr_list.back()->
	      expr_list.front()->int_val ;
	    if(OP_UNARY_MINUS == end->expr_list.back()->op)
	      time_offset = -time_offset ;
	  } else {
	    cerr << "time level should be of form NAME <op> int in expr "
		 << p << endl ;
	    cerr << "got " << end << endl;
	    break ;
	  }
	if(OP_MINUS == end->op)
	  time_offset = -time_offset ;
	
	v.offset = time_offset ;
	
	v.assign = false ;
	if(OP_ASSIGN == end->op) 
	  v.assign = true ;
      }
      break ; 
      
    case OP_FUNC:
      {
	exprList::const_iterator ei ;
	v.name = e->name ;
	for(ei = e->expr_list.begin(); ei != e->expr_list.end(); ++ei) {
	  variable temp = variable(*ei) ;
	  v.v_ids.push_back(temp.ident()) ;
	}
	break ;
      }
    
    case OP_FUNC_BRACE:
      {
	exprList::const_iterator ei ;
	ei = e->expr_list.begin() ;
	exprP func_list = *ei ;
	++ei ;
	exprP brace_list = *ei ;
	string fname = e->name ;
	exprList time_list ;
	if(brace_list->op == OP_COMMA)
	  time_list = brace_list->expr_list ;
	else
	  time_list.push_back(brace_list) ;
	
	exprP time_expr = new expression(OP_NAME_BRACE,fname,time_list) ;
	v.name = fname ;
	v.time_id = time_ident(time_expr) ;
	end = time_expr->expr_list.back() ;
	if(OP_PLUS == end->op  ||
	   OP_MINUS == end->op ||
	   OP_ASSIGN == end->op) {
	  if(end->expr_list.size() != 2) {
	    cerr << "syntax error in time label of name expression "
		 << p << endl ;
	    cerr << "error found while parsing " << end<< endl ;
	    break ;
	  }
	  exprP tname = end->expr_list.front() ;
	  int time_offset = 0 ;
	  if(OP_NAME != tname->op) {
	    cerr << "syntax error in time label of name expression "
		 << p << endl ;
	    cerr << "time level should be of form NAME <op> int" << endl;
	    break ;
	  } else
	    if(OP_INT == end->expr_list.back()->op) {
	      time_offset = end->expr_list.back()->int_val ;
	    } else if((OP_UNARY_PLUS  == end->expr_list.back()->op ||
		       OP_UNARY_MINUS == end->expr_list.back()->op) &&
		      OP_INT == end->expr_list.back()->expr_list.front()->op){
	      time_offset = end->expr_list.back()->
		expr_list.front()->int_val ;
	      if(OP_UNARY_MINUS == end->expr_list.back()->op)
		time_offset = -time_offset ;
	    } else {
	      cerr << "time level should be of form NAME <op> int in expr "
		   << p << endl ;
	      cerr << "got " << end << endl;
	      break ;
	    }
	  if(OP_MINUS == end->op)
	    time_offset = -time_offset ;
	  
	  v.offset = time_offset ;
	  
	  v.assign = false ;
	  if(OP_ASSIGN == end->op) 
	    v.assign = true ;
	}
	
	exprList var_list ;
	if(func_list->op == OP_COMMA)
	  var_list = func_list->expr_list ;
	else
	  var_list.push_back(func_list) ;
	
	for(ei = var_list.begin();ei!=var_list.end();++ei) {
	  variable temp = variable(*ei) ;
	  v.v_ids.push_back(temp.ident()) ;
	}
	
	break ;
      }
    case OP_NIL:
      v.name = e->name ;
      break ;
      
    default:
      cerr << "unable to interpret expression " << e << endl ;
      break ;
    }
    id = vdb->vars.get_id(v) ;
  }
  
  variable::variable(const time_ident &t) {
    info v ;
    v.tvar = true ;
    v.time_id = t ;
    v.name = t.level_name() ;
    v.offset = 0 ;
    v.assign = false ;
    //v.rename = false ;
    create_vdb() ;
    id = vdb->vars.get_id(v) ;
  }
  
  variable::variable(const variable &v, const std::vector<int> &vint) {
    create_vdb() ;
    info v2 = v.get_info() ;
    v2.v_ids = vint ;
    id = vdb->vars.get_id(v2) ;
  }
  
  variable::variable(const variable &v, const time_ident &t) {
    create_vdb() ;
    info v2 = v.get_info() ;
    v2.time_id = t ;
    v2.assign = false ;
    v2.offset = 0 ;
    id = vdb->vars.get_id(v2) ;
  }

  // this constructs a variable with the time_ident prepended
  variable::variable(const time_ident& t, const variable& v) {
    create_vdb() ;
    info v2 = v.get_info() ;
    time_ident new_time = prepend_time(t,v2.time_id) ;
    v2.time_id = new_time ;
    id = vdb->vars.get_id(v2) ;
  }

  variable variable::info::parent() const {
    info vi = *this ;
    vi.time_id = time_id.parent() ;
    return variable(vi) ;
  }
  
  variable variable::info::drop_assign() const {
    info vi = *this ;
    vi.assign = false ;
    return variable(vi) ;
  }
  
  variable variable::info::drop_priority() const {
    info vi = *this ;
    for(int i = 0;i<int(vi.priority.size())-1;++i)
      vi.priority[i] = vi.priority[i+1] ;
    vi.priority.pop_back() ;
    return variable(vi) ;
  }

  variable variable::info::new_offset(int o) const {
    info vi = *this ;
    vi.offset = o ;
    return variable(vi) ;
  }
  
  variable variable::info::drop_namespace() const {
    info vi = *this ;
    if (vi.namespac.empty())
      return variable(vi);
    for(int i = 0;i<int(vi.namespac.size())-1;++i) {
      vi.namespac[i] = vi.namespac[i+1] ;
    }
    vi.namespac.pop_back() ;
    return variable(vi) ;
  }
  
  variable variable::info::add_namespace(const std::string& n) const{
    info vi = *this;
    vi.namespac.insert(vi.namespac.begin(),1,n) ;
    //vi.rename = true ;
    std::vector<int> tmp_vids ;
    for(std::vector<int>::const_iterator vvi = v_ids.begin() ;vvi != v_ids.end(); ++vvi) {
      variable temp = variable(*vvi) ;
      temp = temp.get_info().add_namespace(n) ;
      tmp_vids.push_back(temp.ident()) ;
    }
    vi.v_ids = tmp_vids ;
    return variable(vi);
  }
  /*
  variable variable::info::set_rename() const {
    info vi = *this;
    vi.rename = true ;
    return variable(vi);
  }
  */
  
  variable variable::info::change_time(time_ident ti) const {
    info vi = *this ;
    vi.time_id = ti ;
    return variable(vi) ;
  }
  
ostream &variable::info::Print(ostream &s ) const {
    if(tvar)
      s << "$" ;
    
   
		
    if(priority.begin() != priority.end()) {
      for(vector<string>::const_iterator i =priority.begin();
	  i!=priority.end();++i)
	s<< *i <<"::" ;
      if(namespac.begin() != namespac.end()) {
	for(vector<string>::const_iterator i=namespac.begin();
	    i!=namespac.end();++i) {
	  s << *i << "@" ;
	}
      }
    }
    else {
      if(namespac.begin() != namespac.end()) {
	for(vector<string>::const_iterator i=namespac.begin();
	    i!=namespac.end();++i) {
	  s << *i << "@" ;
	}
      }
    }
    s<< name ;
    
    if(v_ids.size() != 0) {
      s << "(";
      std::vector<int>::const_iterator vi = v_ids.begin() ;
      variable temp = variable(*vi) ;
      temp.get_info().Print(s) ;
      ++vi ;
      for(;vi != v_ids.end(); ++vi) {
	s << ",";
	temp = variable(*vi) ;
	temp.get_info().Print(s) ;
      }
      s << ")" ;
    }
    
    if(time_id != time_ident()) {
      s << "{" << time_id ;
      if(assign)
	s << "=" << offset ;
      else if(offset != 0) {
	if(offset>0) 
	  s << "+" << offset ;
	else
	  s << offset ;
      }
      s << "}" ;
    }
    return s;
  }
  
  variableSet::variableSet(const exprP &e) {
    exprList l = collect_associative_op(e,OP_COMMA) ;
    
    for(exprList::const_iterator i=l.begin();i!=l.end();++i)
      *this += variable(*i) ; 
  }

  vector<string> variableSet::lexico_sort() const {
    vector<string> name ;    
    if(size() == 0)
      return name ;
    for(variableSet::const_iterator vi=begin();vi!=end();++vi) {
      std::ostringstream ss ;
      ss << *vi ;
      name.push_back(ss.str()) ;
    }
    std::sort(name.begin(),name.end(),lexico_cmp) ;
    return name ;
  }

  /* original one, does not support ordered output
  ostream &variableSet::Print(ostream &s) const
  {
    variableSet::const_iterator i = begin() ;
    int sz = size() ;
    if(sz == 0)
      s << "()" ;
    else if(sz == 1)
      s << *i ;
    else {
      s << "(" << *i ;
      ++i ;
      for(;i!=end();++i) 
        s << "," << *i ;
      s << ")" ;
    }
    return s;
  }
  */
  ostream& variableSet::Print(ostream& s) const {
    vector<string> sorted_set = lexico_sort() ;
    vector<string>::const_iterator vi = sorted_set.begin() ;
    int sz = size() ;
    if(sz == 0)
      s << "()" ;
    else if(sz == 1)
      s << *vi ;
    else {
      s << "(" << *vi ;
      ++vi ;
      for(;vi!=sorted_set.end();++vi)
        s << "," << *vi ;
      s << ")" ;
    }
    return s ;
  }

        
  vmap_info::vmap_info(const exprP &e) {
    exprList l = collect_associative_op(e,OP_ARROW) ;
    exprList v = collect_associative_op(l.back(),OP_COMMA) ;

    for(exprList::const_iterator i=v.begin();i!=v.end();++i) {
      
      if((*i)->op == OP_ASSIGN) {
	variable varid((*i)->expr_list.front()) ;
	variable asnid((*i)->expr_list.back()) ;

	var += varid ;
            assign.push_back(make_pair(varid,asnid)) ;
      } else {
	var += variable(*i) ;
      }
    }
    
    l.pop_back() ;
    for(exprList::const_iterator j = l.begin();j!=l.end();++j)
      mapping.push_back(variableSet(*j)) ;
  }

  ostream &vmap_info::Print(ostream &s) const {
    for(size_t j=0;j<mapping.size();++j) 
      s << mapping[j] << "->" ;
    s << var ;
    return s ;
  }
  
}
